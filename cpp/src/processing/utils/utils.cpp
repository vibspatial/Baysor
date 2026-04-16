#include "baysor/processing/utils/utils.h"
#include <third_party/nanoflann.hpp>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <omp.h>

namespace baysor {

// ============================================================================
// count_array_sparse
// ============================================================================

Eigen::SparseVector<float> count_array_sparse(
    const int* values, int n, int total,
    const double* weights,
    bool normalize
) {
    if (n == 0) return Eigen::SparseVector<float>(total);

    // Sort indices by gene ID
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](int a, int b) {
        return values[a] < values[b];
    });

    // Accumulate counts grouped by gene ID
    std::vector<int> indices;
    std::vector<float> counts;
    constexpr double min_val = 1e-5;

    int last_id = values[perm[0]];
    double cnt = 0.0;
    int id = 0;

    for (int pi = 0; pi < n; ++pi) {
        int i = perm[pi];
        id = values[i];
        if (id != last_id) {
            if (cnt > min_val) {
                counts.push_back(static_cast<float>(cnt));
                indices.push_back(last_id - 1); // 1-based gene IDs -> 0-based sparse index
            }
            cnt = 0.0;
            last_id = id;
        }
        cnt += weights ? weights[i] : 1.0;
    }
    if (cnt > min_val) {
        counts.push_back(static_cast<float>(cnt));
        indices.push_back(id - 1); // 1-based -> 0-based
    }

    if (normalize && !counts.empty()) {
        float s = 0.0f;
        for (float c : counts) s += c;
        if (s > 0.0f) {
            for (float& c : counts) c /= s;
        }
    }

    // Build SparseVector
    Eigen::SparseVector<float> sv(total);
    sv.reserve(static_cast<int>(indices.size()));
    for (size_t i = 0; i < indices.size(); ++i) {
        if (indices[i] >= 0 && indices[i] < total) {
            sv.insert(indices[i]) = counts[i];
        }
    }
    return sv;
}

// ============================================================================
// knn_parallel (nanoflann-based)
// ============================================================================

// Adaptor for dims x N Eigen matrix (columns = points)
struct EigenColMajorAdaptor {
    const Eigen::MatrixXd& mat;
    EigenColMajorAdaptor(const Eigen::MatrixXd& m) : mat(m) {}

    inline size_t kdtree_get_point_count() const { return static_cast<size_t>(mat.cols()); }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return mat(static_cast<Eigen::Index>(dim), static_cast<Eigen::Index>(idx));
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, EigenColMajorAdaptor>,
    EigenColMajorAdaptor,
    -1,  // dynamic dimensionality
    int  // index type
>;

KnnResult knn_parallel(
    const Eigen::MatrixXd& tree_points,
    const Eigen::MatrixXd& query_points,
    int k,
    bool sorted
) {
    const int n_dims = static_cast<int>(tree_points.rows());
    const int n_tree = static_cast<int>(tree_points.cols());
    const int n_query = static_cast<int>(query_points.cols());

    if (n_tree == 0 || n_query == 0 || k <= 0) {
        return {};
    }

    // Clamp k to available points
    k = std::min(k, n_tree);

    // Build KD-tree
    EigenColMajorAdaptor adaptor(tree_points);
    KDTree tree(n_dims, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(/* max_leaf = */ 10));

    KnnResult result;
    result.indices.resize(n_query);
    result.distances.resize(n_query);

    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < n_query; ++i) {
        result.indices[i].resize(k);
        result.distances[i].resize(k);

        nanoflann::KNNResultSet<double, int> resultSet(k);
        resultSet.init(result.indices[i].data(), result.distances[i].data());
        tree.findNeighbors(resultSet, query_points.col(i).data());

        // nanoflann returns squared distances — convert to actual distances
        for (int j = 0; j < k; ++j) {
            result.distances[i][j] = std::sqrt(result.distances[i][j]);
        }
    }

    return result;
}

} // namespace baysor
