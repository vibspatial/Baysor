#pragma once

#include "baysor/processing/models/adj_list.h"
#include <Eigen/Dense>
#include <vector>

namespace baysor {

/// Result of MRF-based molecule clustering
struct ClusteringResult {
    Eigen::MatrixXd exprs;              // n_clusters × n_genes expression profiles
    std::vector<int> assignment;         // per-molecule cluster ID (1-based)
    Eigen::MatrixXd assignment_probs;    // n_clusters × n_molecules soft assignments
    std::vector<double> diffs;           // max_diff convergence trace
    std::vector<double> change_fracs;    // fraction-changed convergence trace
};

/// Core MRF-based categorical EM for molecule clustering.
/// exprs_init: optional pre-computed n_clusters × n_genes expression profiles
///             (e.g. from ICA). Pass nullptr to use deterministic hash fallback.
ClusteringResult cluster_molecules_on_mrf(
    const std::vector<int>& genes,      // 1-based gene IDs
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    int n_clusters,
    double tol       = 0.01,
    double mrf_weight = 1.0,
    int max_iters    = -1,              // -1 = auto (max(10000, n_mols/200))
    bool verbose     = true,
    const Eigen::MatrixXd* exprs_init = nullptr  // nullptr → hash-based fallback
);

/// ICA-initialized molecule clustering (matches Julia's DataFrame wrapper).
/// Computes pairwise gene spatial co-occurrence → FastICA → initializes EM.
/// Falls back to hash initialization if ICA fails to converge.
ClusteringResult cluster_molecules_ica(
    const std::vector<int>& genes,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    int n_clusters,
    double tol       = 0.01,
    double mrf_weight = 1.0,
    int max_iters    = -1,
    bool verbose     = true
);

} // namespace baysor
