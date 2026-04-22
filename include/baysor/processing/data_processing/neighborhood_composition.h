#pragma once

#include <vector>
#include <functional>
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

Eigen::SparseMatrix<float> neighborhood_count_matrix_subset(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const std::vector<int>& query_ids,
    int k,
    int n_genes = -1,
    const std::vector<double>* confidences = nullptr,
    bool normalize_by_dist = true,
    bool normalize = true,
    double distance_floor = -1.0
);

double neighborhood_distance_floor(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>* query_ids = nullptr
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

Eigen::MatrixXf project_gene_vectors(
    const Eigen::MatrixXf& gene_emb_t,
    const Eigen::SparseMatrix<float>& count_matrix
);

Eigen::MatrixXf project_neighborhood_vectors(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int k,
    const Eigen::MatrixXf& gene_emb_t,
    int n_genes = -1,
    const std::vector<int>* query_ids = nullptr,
    const std::vector<double>* confidences = nullptr,
    bool normalize_by_dist = true,
    bool normalize = true,
    double distance_floor = -1.0,
    bool log_transform = false
);

void stream_projected_neighborhood_vectors(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int k,
    const Eigen::MatrixXf& gene_emb_t,
    int n_genes = -1,
    const std::vector<int>* query_ids = nullptr,
    const std::vector<double>* confidences = nullptr,
    bool normalize_by_dist = true,
    bool normalize = true,
    double distance_floor = -1.0,
    bool log_transform = false,
    int block_size = 32768,
    const std::function<void(int, const std::vector<int>&, const Eigen::MatrixXf&)>& callback = nullptr
);

} // namespace baysor
