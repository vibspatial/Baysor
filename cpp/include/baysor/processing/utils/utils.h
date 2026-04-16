#pragma once

#include <Eigen/Sparse>
#include <vector>

namespace baysor {

/// Build a sparse count vector from gene IDs (like Julia's count_array_sparse)
Eigen::SparseVector<float> count_array_sparse(
    const int* values, int n, int total,
    const double* weights = nullptr,
    bool normalize = false
);

/// Parallel KNN query: for each column of query_points, find k nearest neighbors in tree_points.
/// Returns (indices, distances), each n_points x k.
struct KnnResult {
    std::vector<std::vector<int>> indices;       // [n_points][k]
    std::vector<std::vector<double>> distances;  // [n_points][k]
};

KnnResult knn_parallel(
    const Eigen::MatrixXd& tree_points,
    const Eigen::MatrixXd& query_points,
    int k,
    bool sorted = false
);

} // namespace baysor
