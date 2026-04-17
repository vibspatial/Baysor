#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace baysor {

/// Build sparse neighborhood count matrix: n_genes x n_molecules
/// For each molecule, counts genes among its K nearest neighbors.
Eigen::SparseMatrix<float> neighborhood_count_matrix(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int k,
    int n_genes = -1,
    const std::vector<double>* confidences = nullptr,
    bool normalize_by_dist = true,
    bool normalize = true
);

/// Estimate gene vectors via randomized indexing.
/// Returns n_components x n_molecules (per_molecule=true) or n_components x n_genes matrix.
/// Kept as float throughout to avoid materializing a dense copy of the large sparse input.
Eigen::MatrixXf estimate_gene_vectors(
    const Eigen::SparseMatrix<float>& count_matrix,
    const std::vector<int>& gene_ids,
    int n_components,
    const std::string& method = "ri",
    bool per_molecule = false
);

} // namespace baysor
