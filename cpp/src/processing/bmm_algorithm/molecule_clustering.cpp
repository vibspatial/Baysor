#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/utils/general.h"

#include <Eigen/Dense>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace baysor {

// ============================================================================
// fast_ica — symmetric deflation FastICA (tanh nonlinearity)
//
// Input:  X  — data matrix (n_features × n_samples)
//         n_components — number of independent components
// Returns: mixing matrix A (n_features × n_components)
//          where X ≈ A * S  (columns of A are the mixing directions)
//
// Mirrors Julia's MultivariateStats.ICA usage in cluster_molecules_on_mrf.
// ============================================================================
static Eigen::MatrixXd fast_ica(
    const Eigen::MatrixXd& X,
    int n_components,
    int max_iter = 1000,
    double tol   = 1e-5,
    unsigned int seed = 42
) {
    int n_features = static_cast<int>(X.rows());
    int n_samples  = static_cast<int>(X.cols());
    if (n_components > n_features) n_components = n_features;

    // 1. Center: subtract column mean
    Eigen::VectorXd mean_vec = X.rowwise().mean();
    Eigen::MatrixXd Xc = X.colwise() - mean_vec;

    // 2. Whiten: PCA → keep top n_components directions
    Eigen::MatrixXd cov = (Xc * Xc.transpose()) / static_cast<double>(n_samples);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
    // eigenvalues in ascending order; take last n_components
    Eigen::VectorXd lambdas = eig.eigenvalues().tail(n_components);
    Eigen::MatrixXd V       = eig.eigenvectors().rightCols(n_components);  // n_features × n_components

    // Clamp tiny eigenvalues to avoid division by ~0
    for (int i = 0; i < n_components; ++i)
        if (lambdas(i) < 1e-10) lambdas(i) = 1e-10;

    // Whitened data Z: n_components × n_samples
    // Z = diag(1/sqrt(λ)) * V^T * Xc
    Eigen::MatrixXd Z = (lambdas.cwiseSqrt().cwiseInverse().asDiagonal()) * (V.transpose() * Xc);

    // 3. Fixed-point FastICA (deflation)
    // Each row of W_z is an unmixing direction in the whitened space
    Eigen::MatrixXd W_z(n_components, n_components);
    W_z.setZero();

    std::mt19937 rng(seed);
    std::normal_distribution<double> ndist(0.0, 1.0);

    for (int i = 0; i < n_components; ++i) {
        // Random init
        Eigen::VectorXd w(n_components);
        for (int j = 0; j < n_components; ++j) w(j) = ndist(rng);
        w.normalize();

        for (int iter = 0; iter < max_iter; ++iter) {
            Eigen::VectorXd proj = Z.transpose() * w;           // n_samples
            Eigen::VectorXd g    = proj.array().tanh();          // g(u)  = tanh(u)
            Eigen::VectorXd gd   = 1.0 - g.array().square();    // g'(u) = 1 - tanh²(u)

            Eigen::VectorXd w_new = (Z * g) / n_samples
                                  - gd.mean() * w;

            // Gram-Schmidt deflation
            for (int k = 0; k < i; ++k)
                w_new -= w_new.dot(W_z.row(k)) * W_z.row(k).transpose();

            w_new.normalize();

            double delta = (w_new - w).norm();
            w = w_new;
            if (delta < tol) break;
        }

        W_z.row(i) = w.transpose();
    }

    // 4. Compute mixing matrix in original space:
    //    A_orig = V * diag(sqrt(λ)) * W_z^T
    Eigen::MatrixXd A = V * lambdas.cwiseSqrt().asDiagonal() * W_z.transpose();
    // A: n_features × n_components
    return A;
}

// ============================================================================
// cluster_molecules_on_mrf — categorical variant
// ============================================================================
//
// Port of Julia's cluster_molecules_on_mrf (molecule_clustering.jl)
// using the CatMixture (categorical gene-expression) mixture model.
//
// Layout conventions (matching Julia):
//   exprs            : n_clusters × n_genes    (row = cluster, col = gene)
//   assignment_probs : n_clusters × n_molecules (row = cluster, col = mol)
//   genes            : 1-based (0 = missing/unknown gene, skip)
//
// Optional exprs_init: pre-computed n_clusters × n_genes expression profiles
// (e.g. from ICA). When nullptr, falls back to deterministic hash perturbation.

ClusteringResult cluster_molecules_on_mrf(
    const std::vector<int>& genes,   // 1-based gene IDs
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    int n_clusters,
    double tol,
    double mrf_weight,
    int max_iters,
    bool verbose,
    const Eigen::MatrixXd* exprs_init
) {
    int n_mols  = static_cast<int>(genes.size());
    if (n_mols == 0 || n_clusters <= 1) return {};

    // Infer n_genes from max gene ID (genes are 1-based)
    int n_genes = 0;
    for (int g : genes) if (g > n_genes) n_genes = g;
    if (n_genes == 0) return {};

    if (max_iters <= 0)
        max_iters = std::max(10000, n_mols / 200);

    constexpr int n_iters_without_update = 20;

    // ------------------------------------------------------------------
    // Step 0: Pre-multiply adj weights by neighbor confidence
    //   adj_weights_conf[i][j] = original_weight[i][j] * confidence[neighbor[j]]
    // Stored as a flat mirror of adj_list (same indptr).
    // ------------------------------------------------------------------
    std::vector<double> adj_w_conf(adj_list.nnz());
    for (int i = 0; i < n_mols; ++i) {
        int start = adj_list.indptr[i];
        int end   = adj_list.indptr[i + 1];
        const int32_t* nb  = adj_list.indices.data() + start;
        const double*  wt  = adj_list.weights.data()  + start;
        double*        out = adj_w_conf.data()         + start;
        for (int j = 0; j < end - start; ++j) {
            out[j] = wt[j] * confidence[nb[j]];
        }
    }

    // ------------------------------------------------------------------
    // Step 1: Initialize expression profiles
    //
    // If exprs_init is provided (e.g. from ICA), use it directly (after
    // normalizing rows to sum to 1).  Otherwise, fall back to Julia's
    // deterministic hash-perturbed gene-frequency initialization
    // (init_cell_type_exprs with init_mod=10000).
    // ------------------------------------------------------------------

    // Global gene frequency (0-based genes) — needed for fallback init
    std::vector<double> gene_freq(n_genes, 0.0);
    int n_valid = 0;
    for (int g1b : genes) {
        int g0 = g1b - 1;
        if (g0 >= 0 && g0 < n_genes) { gene_freq[g0] += 1.0; ++n_valid; }
    }
    if (n_valid > 0) for (double& f : gene_freq) f /= n_valid;

    Eigen::MatrixXd exprs(n_clusters, n_genes);

    if (exprs_init && exprs_init->rows() == n_clusters
                   && exprs_init->cols() == n_genes) {
        // Use provided initialization (ICA result), row-normalize
        exprs = *exprs_init;
        for (int k = 0; k < n_clusters; ++k) {
            double rs = exprs.row(k).sum();
            if (rs > 1e-10) exprs.row(k) /= rs;
            else            exprs.row(k).fill(1.0 / n_genes);
        }
    } else {
        // Fallback: deterministic hash perturbation of global gene frequency
        // Mirrors Julia's init_cell_type_exprs (init_mod=10000 path)
        constexpr int init_mod = 10000;
        for (int k = 0; k < n_clusters; ++k) {
            double row_sum = 0.0;
            for (int g = 0; g < n_genes; ++g) {
                std::size_t h = std::hash<long long>{}(
                    static_cast<long long>(g + 1) *
                    static_cast<long long>((k + 1) * (k + 1)));
                double noise = static_cast<double>(h % init_mod) / 100000.0;
                exprs(k, g) = gene_freq[g] * (0.95 + noise);
                row_sum += exprs(k, g);
            }
            // Pseudocount normalize matching Julia: (x+1)/(sum+1) per element
            double norm = row_sum + static_cast<double>(n_genes);
            for (int g = 0; g < n_genes; ++g)
                exprs(k, g) = (exprs(k, g) + 1.0) / (norm + n_genes);
            double rs2 = exprs.row(k).sum();
            if (rs2 > 0) exprs.row(k) /= rs2;
        }
    }

    // assignment_probs[k][i] = exprs[k][gene[i]], then column-normalize
    // (matches Julia's init_assignment_probs_inner)
    Eigen::MatrixXd probs(n_clusters, n_mols);
    for (int i = 0; i < n_mols; ++i) {
        int g0 = genes[i] - 1;
        double col_sum = 0.0;
        for (int k = 0; k < n_clusters; ++k) {
            double p = (g0 >= 0) ? exprs(k, g0) : (1.0 / n_clusters);
            probs(k, i) = p;
            col_sum += p;
        }
        if (col_sum > 1e-100)
            for (int k = 0; k < n_clusters; ++k) probs(k, i) /= col_sum;
        else
            for (int k = 0; k < n_clusters; ++k) probs(k, i) = 1.0 / n_clusters;
    }

    Eigen::MatrixXd prev_probs(n_clusters, n_mols);

    std::vector<double> max_diffs;
    std::vector<double> change_fracs;
    max_diffs.reserve(max_iters);
    change_fracs.reserve(max_iters);

    // ------------------------------------------------------------------
    // Step 2: EM loop
    // ------------------------------------------------------------------
    int n_iters_done = 0;
    for (int iter = 0; iter < max_iters; ++iter) {
        n_iters_done = iter + 1;
        prev_probs = probs;

        // ---- E-step (parallel over molecules) ----
        #pragma omp parallel for schedule(dynamic, 512)
        for (int i = 0; i < n_mols; ++i) {
            int  g0    = genes[i] - 1;  // 0-based gene (< 0 if missing)
            int  start = adj_list.indptr[i];
            int  end   = adj_list.indptr[i + 1];
            const int32_t* nb_ids = adj_list.indices.data() + start;
            const double*  nb_wt  = adj_w_conf.data()        + start;
            int  n_nb  = end - start;

            double col_sum = 0.0;
            for (int k = 0; k < n_clusters; ++k) {
                // MRF term: weighted sum of neighbor probabilities for cluster k
                double c_d = 0.0;
                for (int j = 0; j < n_nb; ++j) {
                    double a_p = prev_probs(k, nb_ids[j]);
                    if (a_p > 1e-5) c_d += nb_wt[j] * a_p;
                }
                double mrf_prior = std::exp(mrf_weight * c_d);

                // Expression likelihood (skip if gene unknown)
                double expr_ll = (g0 >= 0) ? exprs(k, g0) : 1.0;
                probs(k, i) = expr_ll * mrf_prior;
                col_sum    += probs(k, i);
            }

            // Normalize column
            if (col_sum > 1e-100) {
                for (int k = 0; k < n_clusters; ++k) probs(k, i) /= col_sum;
            } else {
                for (int k = 0; k < n_clusters; ++k) probs(k, i) = 1.0 / n_clusters;
            }
        }

        // ---- M-step with pseudocount ----
        exprs.setZero();
        for (int i = 0; i < n_mols; ++i) {
            int g0 = genes[i] - 1;
            if (g0 < 0) continue;
            double conf = confidence[i];
            for (int k = 0; k < n_clusters; ++k) {
                exprs(k, g0) += conf * probs(k, i);
            }
        }
        // Add pseudocount and row-normalize
        exprs.array() += 1.0;
        for (int k = 0; k < n_clusters; ++k) {
            double row_sum = exprs.row(k).sum();
            if (row_sum > 0) exprs.row(k) /= row_sum;
        }

        // ---- Convergence check ----
        double max_diff = 0.0;
        int n_changed = 0;
        for (int i = 0; i < n_mols; ++i) {
            double conf = confidence[i];
            double mol_max = 0.0;
            for (int k = 0; k < n_clusters; ++k) {
                double d = std::abs(probs(k, i) - prev_probs(k, i)) * conf;
                if (d > mol_max) mol_max = d;
                if (d > max_diff) max_diff = d;
            }
            if (mol_max > 0.05) ++n_changed;
        }
        max_diffs.push_back(max_diff);
        change_fracs.push_back(static_cast<double>(n_changed) / n_mols);

        if (verbose && (iter % 100 == 0 || iter < 5)) {
            spdlog::info("  Clustering iter {:4d}: max_diff={:.4f}, change_frac={:.4f}",
                         iter + 1, max_diff, change_fracs.back());
        }

        // Stop if last n_iters_without_update all below tol
        if (iter >= n_iters_without_update) {
            int look_back = std::min(n_iters_without_update,
                                     static_cast<int>(max_diffs.size()));
            double worst = 0.0;
            for (int t = static_cast<int>(max_diffs.size()) - look_back;
                 t < static_cast<int>(max_diffs.size()); ++t) {
                if (max_diffs[t] > worst) worst = max_diffs[t];
            }
            if (worst < tol) {
                if (verbose)
                    spdlog::info("Clustering converged after {} iterations. Max diff: {:.4f}",
                                 iter + 1, max_diff);
                break;
            }
        }
    }

    if (verbose && !max_diffs.empty()) {
        spdlog::info("Clustering stopped after {} iterations. Max diff: {:.4f}. Converged: {}",
                     n_iters_done, max_diffs.back(), max_diffs.back() < tol ? "true" : "false");
    }

    // ---- Final M-step without pseudocount ----
    exprs.setZero();
    for (int i = 0; i < n_mols; ++i) {
        int g0 = genes[i] - 1;
        if (g0 < 0) continue;
        double conf = confidence[i];
        for (int k = 0; k < n_clusters; ++k) {
            exprs(k, g0) += conf * probs(k, i);
        }
    }
    for (int k = 0; k < n_clusters; ++k) {
        double row_sum = exprs.row(k).sum();
        if (row_sum > 0) exprs.row(k) /= row_sum;
    }

    // ---- Hard assignment: 1-based cluster IDs ----
    std::vector<int> assignment(n_mols);
    for (int i = 0; i < n_mols; ++i) {
        int best_k = 0;
        double best_p = probs(0, i);
        for (int k = 1; k < n_clusters; ++k) {
            if (probs(k, i) > best_p) { best_p = probs(k, i); best_k = k; }
        }
        assignment[i] = best_k + 1;  // 1-based
    }

    return ClusteringResult{
        std::move(exprs),
        std::move(assignment),
        std::move(probs),
        std::move(max_diffs),
        std::move(change_fracs)
    };
}

// ============================================================================
// cluster_molecules_ica
//
// Port of Julia's DataFrame-based wrapper for cluster_molecules_on_mrf:
//   1. Compute pairwise gene spatial co-occurrence (already in color_utils)
//   2. Run FastICA on the correlation matrix to get n_clusters gene profiles
//   3. Call core EM with ICA init; fall back to hash init if ICA throws
// ============================================================================
ClusteringResult cluster_molecules_ica(
    const std::vector<int>& genes,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    int n_clusters,
    double tol,
    double mrf_weight,
    int max_iters,
    bool verbose
) {
    if (n_clusters <= 1) return {};

    // Infer n_genes
    int n_genes = 0;
    for (int g : genes) if (g > n_genes) n_genes = g;

    // 1. Gene spatial correlation matrix (n_genes × n_genes)
    Eigen::MatrixXd cor_mat = pairwise_gene_spatial_cor(genes, confidence, adj_list);
    // cor_mat is indexed 0-based (genes 0..n_genes-1 from 1-based input)

    // 2. FastICA on the correlation matrix → mixing matrix (n_genes × n_clusters)
    std::unique_ptr<Eigen::MatrixXd> exprs_init_ptr;
    try {
        Eigen::MatrixXd A = fast_ica(cor_mat, n_clusters);
        // Convert mixing matrix to expression profiles:
        //   ct_exprs_init[k][g] = abs(A[g][k]) / sum_g'(abs(A[g'][k]))
        // Matches Julia: (abs.(ica_fit.W) ./ sum(abs.(ica_fit.W), dims=1))'
        Eigen::MatrixXd exprs(n_clusters, n_genes);
        for (int k = 0; k < n_clusters; ++k) {
            double col_sum = A.col(k).cwiseAbs().sum();
            if (col_sum < 1e-10) col_sum = 1.0;
            for (int g = 0; g < n_genes; ++g)
                exprs(k, g) = std::abs(A(g, k)) / col_sum;
        }
        exprs_init_ptr = std::make_unique<Eigen::MatrixXd>(std::move(exprs));
        if (verbose) spdlog::info("ICA initialization succeeded ({} components).", n_clusters);
    } catch (const std::exception& e) {
        spdlog::warn("ICA did not converge ({}), falling back to hash initialization.", e.what());
    } catch (...) {
        spdlog::warn("ICA failed, falling back to hash initialization.");
    }

    // 3. Core EM with ICA init (or nullptr → hash fallback)
    return cluster_molecules_on_mrf(
        genes, adj_list, confidence,
        n_clusters, tol, mrf_weight, max_iters, verbose,
        exprs_init_ptr.get()
    );
}

} // namespace baysor
