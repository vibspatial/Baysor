#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/processing/utils/utils.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <omp.h>

namespace baysor {

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

    if (n_genes == 0 || n_mols == 0) return Eigen::MatrixXf();

    // Random projection matrix: n_genes x n_components (float, matches count_matrix)
    std::mt19937 rng(42);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    Eigen::MatrixXf random_vectors(n_genes, n_components);
    for (int i = 0; i < n_genes; ++i)
        for (int j = 0; j < n_components; ++j)
            random_vectors(i, j) = normal(rng);

    // Co-expression matrix: sparse(84xN) * sparse(84xN)^T → sparse(84x84)
    // Materialize to dense immediately since it's small (n_genes x n_genes).
    // count_matrix stays SparseMatrix<float> — no dense copy is made.
    Eigen::SparseMatrix<float> coexpr_sp = count_matrix * count_matrix.transpose();
    Eigen::MatrixXf coexpr(coexpr_sp);  // 84×84 dense, cheap

    // Variance clipping (var_clip=0.05): cap diagonal so highly self-correlated
    // genes don't dominate (matches Julia's generate_randomized_gene_vectors).
    constexpr float var_clip = 0.05f;
    if (var_clip > 0 && n_genes > 1) {
        Eigen::VectorXf diag_vals = coexpr.diagonal();
        Eigen::VectorXf total_var = coexpr.rowwise().sum();

        std::vector<float> df_sorted(n_genes);
        for (int i = 0; i < n_genes; ++i)
            df_sorted[i] = (total_var(i) > 0) ? diag_vals(i) / total_var(i) : 0.0f;
        std::sort(df_sorted.begin(), df_sorted.end());
        float q = df_sorted[static_cast<int>((1.0f - var_clip) * (n_genes - 1))];

        for (int i = 0; i < n_genes; ++i)
            coexpr(i, i) = std::min(q * total_var(i), diag_vals(i));
    }

    // Gene embedding: (coexpr * random_vectors) / row_sums — all 84×20, trivial
    Eigen::VectorXf row_sums = coexpr.rowwise().sum();
    Eigen::MatrixXf gene_emb = coexpr * random_vectors;
    for (int i = 0; i < n_genes; ++i)
        if (row_sums(i) > 0) gene_emb.row(i) /= row_sums(i);

    if (per_molecule) {
        // sparse(84xN)^T * dense(84x20) → dense(Nx20), then transpose to 20xN.
        // Eigen dispatches to sparse_time_dense_product, iterating only nonzeros.
        // count_matrix is never densified.
        return (count_matrix.transpose() * gene_emb).transpose();
    }

    return gene_emb.transpose();  // n_components x n_genes
}

} // namespace baysor
