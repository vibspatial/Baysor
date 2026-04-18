#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/processing/utils/utils.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <omp.h>

namespace baysor {

namespace {

Eigen::MatrixXf dense_times_sparse(
    const Eigen::MatrixXf& left,
    const Eigen::SparseMatrix<float>& right
) {
    const int n_components = static_cast<int>(left.rows());
    const int n_cols = static_cast<int>(right.cols());
    Eigen::MatrixXf out = Eigen::MatrixXf::Zero(n_components, n_cols);

    #pragma omp parallel for schedule(dynamic, 256)
    for (int col = 0; col < n_cols; ++col) {
        auto out_col = out.col(col);
        for (Eigen::SparseMatrix<float>::InnerIterator it(right, col); it; ++it) {
            out_col.noalias() += it.value() * left.col(it.row());
        }
    }
    return out;
}

Eigen::MatrixXf dense_times_sparse_transpose(
    const Eigen::MatrixXf& left,
    const Eigen::SparseMatrix<float>& right
) {
    const int n_components = static_cast<int>(left.rows());
    const int n_rows = static_cast<int>(right.rows());
    const int n_cols = static_cast<int>(right.cols());

    const int n_threads = omp_get_max_threads();
    std::vector<Eigen::MatrixXf> locals(
        n_threads, Eigen::MatrixXf::Zero(n_components, n_rows)
    );

    #pragma omp parallel for schedule(dynamic, 256)
    for (int col = 0; col < n_cols; ++col) {
        int tid = omp_get_thread_num();
        auto left_col = left.col(col);
        auto& local = locals[tid];
        for (Eigen::SparseMatrix<float>::InnerIterator it(right, col); it; ++it) {
            local.col(it.row()).noalias() += it.value() * left_col;
        }
    }

    Eigen::MatrixXf out = Eigen::MatrixXf::Zero(n_components, n_rows);
    for (const auto& local : locals) out += local;
    return out;
}

void compute_diag_and_total_var(
    const Eigen::SparseMatrix<float>& count_matrix,
    Eigen::VectorXf& diag_vals,
    Eigen::VectorXf& total_var
) {
    const int n_genes = static_cast<int>(count_matrix.rows());
    const int n_mols = static_cast<int>(count_matrix.cols());

    std::vector<float> col_sums(n_mols, 0.0f);
    #pragma omp parallel for schedule(dynamic, 256)
    for (int col = 0; col < n_mols; ++col) {
        float sum = 0.0f;
        for (Eigen::SparseMatrix<float>::InnerIterator it(count_matrix, col); it; ++it) {
            sum += it.value();
        }
        col_sums[col] = sum;
    }

    const int n_threads = omp_get_max_threads();
    std::vector<Eigen::VectorXf> diag_locals(n_threads, Eigen::VectorXf::Zero(n_genes));
    std::vector<Eigen::VectorXf> total_locals(n_threads, Eigen::VectorXf::Zero(n_genes));

    #pragma omp parallel for schedule(dynamic, 256)
    for (int col = 0; col < n_mols; ++col) {
        int tid = omp_get_thread_num();
        auto& diag_local = diag_locals[tid];
        auto& total_local = total_locals[tid];
        const float col_sum = col_sums[col];
        for (Eigen::SparseMatrix<float>::InnerIterator it(count_matrix, col); it; ++it) {
            const int row = it.row();
            const float value = it.value();
            diag_local(row) += value * value;
            total_local(row) += value * col_sum;
        }
    }

    diag_vals = Eigen::VectorXf::Zero(n_genes);
    total_var = Eigen::VectorXf::Zero(n_genes);
    for (int t = 0; t < n_threads; ++t) {
        diag_vals += diag_locals[t];
        total_var += total_locals[t];
    }
}

} // namespace

// ============================================================================
// neighborhood_count_matrix
// ============================================================================

Eigen::SparseMatrix<float> neighborhood_count_matrix(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int k,
    int n_genes,
    const std::vector<double>* confidences,
    bool normalize_by_dist,
    bool normalize
) {
    int n = static_cast<int>(pos_data.cols());
    if (n == 0) return Eigen::SparseMatrix<float>(n_genes, 0);

    if (k < 3) k = 3;
    if (!normalize) normalize_by_dist = false;
    if (confidences && normalize_by_dist) {
        confidences = nullptr; // Not supported together
    }

    k = std::min(k, n);
    if (n_genes <= 0) {
        n_genes = *std::max_element(genes.begin(), genes.end());
    }

    // KNN search
    auto knn = knn_parallel(pos_data, pos_data, k, true);

    // For normalize_by_dist: median of closest non-zero distance
    double med_closest_dist = 1e-15;
    if (normalize_by_dist) {
        std::vector<double> closest;
        closest.reserve(n);
        for (int i = 0; i < n; ++i) {
            for (double d : knn.distances[i]) {
                if (d > 1e-15) {
                    closest.push_back(d);
                    break;
                }
            }
        }
        if (!closest.empty()) {
            std::sort(closest.begin(), closest.end());
            med_closest_dist = closest[closest.size() / 2];
        }
    }

    // Build sparse column vectors in parallel, then assemble
    std::vector<Eigen::SparseVector<float>> cols(n);

    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < n; ++i) {
        int ki = static_cast<int>(knn.indices[i].size());
        std::vector<int> neighbor_genes(ki);
        for (int j = 0; j < ki; ++j) {
            neighbor_genes[j] = genes[knn.indices[i][j]];
        }

        if (normalize_by_dist) {
            std::vector<double> weights(ki);
            for (int j = 0; j < ki; ++j) {
                weights[j] = 1.0 / std::max(knn.distances[i][j], med_closest_dist);
            }
            cols[i] = count_array_sparse(neighbor_genes.data(), ki, n_genes, weights.data(), normalize);
        } else if (confidences) {
            std::vector<double> weights(ki);
            for (int j = 0; j < ki; ++j) {
                weights[j] = (*confidences)[knn.indices[i][j]];
            }
            cols[i] = count_array_sparse(neighbor_genes.data(), ki, n_genes, weights.data(), normalize);
        } else {
            cols[i] = count_array_sparse(neighbor_genes.data(), ki, n_genes, nullptr, normalize);
        }
    }

    // Assemble into sparse matrix (n_genes x n_molecules)
    // Count total non-zeros
    int total_nnz = 0;
    for (int i = 0; i < n; ++i) total_nnz += cols[i].nonZeros();

    Eigen::SparseMatrix<float> result(n_genes, n);
    result.reserve(total_nnz);
    for (int i = 0; i < n; ++i) {
        for (Eigen::SparseVector<float>::InnerIterator it(cols[i]); it; ++it) {
            result.insert(it.index(), i) = it.value();
        }
    }
    result.makeCompressed();
    return result;
}

// ============================================================================
// estimate_gene_vectors (randomized indexing)
// ============================================================================

Eigen::MatrixXf estimate_gene_vectors(
    const Eigen::SparseMatrix<float>& count_matrix,
    const std::vector<int>& gene_ids,
    int n_components,
    const std::string& method,
    bool per_molecule
) {
    int n_genes = static_cast<int>(count_matrix.rows());
    int n_mols = static_cast<int>(count_matrix.cols());
    (void)gene_ids;
    (void)method;

    if (n_genes == 0 || n_mols == 0) return Eigen::MatrixXf();

    // Random projection matrix, stored as components x genes so dense_times_sparse
    // can work column-wise against the sparse count matrix.
    std::mt19937 rng(42);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    Eigen::MatrixXf random_vectors_t(n_components, n_genes);
    for (int i = 0; i < n_genes; ++i)
        for (int j = 0; j < n_components; ++j)
            random_vectors_t(j, i) = normal(rng);

    // Avoid explicitly forming coexpr = count_matrix * count_matrix.transpose().
    // We only need:
    //   1. diag(coexpr)
    //   2. row_sums(coexpr)
    //   3. coexpr * random_vectors
    // All of those can be derived by streaming over sparse columns.
    Eigen::VectorXf diag_vals;
    Eigen::VectorXf total_var;
    compute_diag_and_total_var(count_matrix, diag_vals, total_var);

    constexpr float var_clip = 0.05f;
    Eigen::VectorXf clipped_diag = diag_vals;
    if (var_clip > 0 && n_genes > 1) {
        std::vector<float> df_sorted(n_genes);
        for (int i = 0; i < n_genes; ++i)
            df_sorted[i] = (total_var(i) > 0) ? diag_vals(i) / total_var(i) : 0.0f;
        std::sort(df_sorted.begin(), df_sorted.end());
        float q = df_sorted[static_cast<int>((1.0f - var_clip) * (n_genes - 1))];

        for (int i = 0; i < n_genes; ++i)
            clipped_diag(i) = std::min(q * total_var(i), diag_vals(i));
    }

    Eigen::VectorXf delta = diag_vals - clipped_diag;
    Eigen::VectorXf row_sums = total_var - delta;

    // Unclipped coexpr * random_vectors = count_matrix * (count_matrix^T * random_vectors).
    Eigen::MatrixXf proj_mol = dense_times_sparse(random_vectors_t, count_matrix); // k x n_mols
    Eigen::MatrixXf gene_emb_t = dense_times_sparse_transpose(proj_mol, count_matrix); // k x n_genes

    // Diagonal clipping only changes the diagonal of coexpr, so subtract the removed
    // diagonal mass directly from each gene vector, then normalize by the clipped row sums.
    for (int i = 0; i < n_genes; ++i) {
        gene_emb_t.col(i).noalias() -= delta(i) * random_vectors_t.col(i);
        if (row_sums(i) > 0) gene_emb_t.col(i) /= row_sums(i);
    }

    if (per_molecule) {
        return dense_times_sparse(gene_emb_t, count_matrix);
    }

    return gene_emb_t;  // n_components x n_genes
}

} // namespace baysor
