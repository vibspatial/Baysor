#include "baysor/processing/data_processing/triangulation.h"
#include "baysor/processing/utils/utils.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <set>

namespace baysor {

// ============================================================================
// normalize_points
// ============================================================================

Eigen::MatrixXd normalize_points(const Eigen::MatrixXd& points) {
    const int dims = static_cast<int>(points.rows());
    const int n = static_cast<int>(points.cols());
    if (n == 0) return points;

    Eigen::MatrixXd out = points;

    // Shift: subtract per-row minimum
    for (int d = 0; d < dims; ++d) {
        double mn = out.row(d).minCoeff();
        out.row(d).array() -= mn;
    }

    // Scale: divide by global max * 1.1, then shift to [1.01, ~1.9]
    double mx = out.maxCoeff();
    if (mx > 0) {
        out /= mx * 1.1;
    }
    out.array() += 1.01;

    // Jitter duplicate points (KNN with distance=0 creates self-loops)
    if (n > 1) {
        auto knn = knn_parallel(out, out, 2, true);
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> jitter(-1e-5, 1e-5);
        for (int i = 0; i < n; ++i) {
            if (knn.distances[i].size() >= 2 && knn.distances[i][1] < 1e-6) {
                for (int d = 0; d < dims; ++d) {
                    out(d, i) += jitter(rng);
                }
            }
        }
    }

    return out;
}

// ============================================================================
// filter_long_edges
// ============================================================================

void filter_long_edges(AdjacencyResult& result, double n_mads) {
    int n = static_cast<int>(result.edge_dists.size());
    if (n == 0) return;

    // Compute log10 distances
    std::vector<double> log_dists(n);
    for (int i = 0; i < n; ++i) {
        log_dists[i] = std::log10(std::max(result.edge_dists[i], 1e-30));
    }

    // Median
    std::vector<double> sorted_log = log_dists;
    std::sort(sorted_log.begin(), sorted_log.end());
    double median;
    if (n % 2 == 0) {
        median = (sorted_log[n / 2 - 1] + sorted_log[n / 2]) / 2.0;
    } else {
        median = sorted_log[n / 2];
    }

    // MAD (median absolute deviation), normalized
    std::vector<double> abs_devs(n);
    for (int i = 0; i < n; ++i) {
        abs_devs[i] = std::abs(log_dists[i] - median);
    }
    std::sort(abs_devs.begin(), abs_devs.end());
    double mad;
    if (n % 2 == 0) {
        mad = (abs_devs[n / 2 - 1] + abs_devs[n / 2]) / 2.0;
    } else {
        mad = abs_devs[n / 2];
    }
    mad *= 1.4826; // normalize=true (consistent estimator of std)

    double threshold = median + n_mads * mad;

    // Filter
    AdjacencyResult filtered;
    for (int i = 0; i < n; ++i) {
        if (log_dists[i] < threshold) {
            filtered.edge_src.push_back(result.edge_src[i]);
            filtered.edge_dst.push_back(result.edge_dst[i]);
            filtered.edge_dists.push_back(result.edge_dists[i]);
        }
    }
    result = std::move(filtered);
}

// ============================================================================
// adjacency_list
// ============================================================================

AdjacencyResult adjacency_list(
    const Eigen::MatrixXd& points,
    bool filter,
    double n_mads,
    int k_adj,
    AdjacencyType type
) {
    const int dims = static_cast<int>(points.rows());
    const int n = static_cast<int>(points.cols());

    if (n <= 1) return {};

    // Auto-select type
    if (type == AdjacencyType::Auto) {
        type = (dims == 3) ? AdjacencyType::Knn : AdjacencyType::Triangulation;
    }
    if (dims == 3 && type == AdjacencyType::Triangulation) {
        type = AdjacencyType::Knn; // 3D only supports KNN
    }

    Eigen::MatrixXd norm_pts = normalize_points(points);

    // Use a set of ordered pairs to deduplicate edges
    std::set<std::pair<int, int>> edge_set;

    // --- Triangulation (2D only, using CGAL with vertex info) ---
    if (dims == 2 && (type == AdjacencyType::Triangulation || type == AdjacencyType::Both)) {
        // We'll use CGAL with vertex info to store our original index
        using K = CGAL::Exact_predicates_inexact_constructions_kernel;

        // Use Triangulation_vertex_base_with_info to store the original index
        using Vb = CGAL::Triangulation_vertex_base_with_info_2<int, K>;
        using Fb = CGAL::Triangulation_face_base_2<K>;
        using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
        using DT = CGAL::Delaunay_triangulation_2<K, Tds>;
        using Point = K::Point_2;

        std::vector<std::pair<Point, int>> indexed_pts(n);
        for (int i = 0; i < n; ++i) {
            indexed_pts[i] = {Point(norm_pts(0, i), norm_pts(1, i)), i};
        }

        DT dt;
        dt.insert(indexed_pts.begin(), indexed_pts.end());

        // Extract edges
        for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
            auto face = eit->first;
            int idx = eit->second;
            auto v1 = face->vertex((idx + 1) % 3);
            auto v2 = face->vertex((idx + 2) % 3);
            int i1 = v1->info();
            int i2 = v2->info();
            int lo = std::min(i1, i2);
            int hi = std::max(i1, i2);
            edge_set.insert({lo, hi});
        }
    }

    // --- KNN ---
    if (type == AdjacencyType::Knn || type == AdjacencyType::Both) {
        auto knn = knn_parallel(norm_pts, norm_pts, k_adj + 1, true);
        for (int i = 0; i < n; ++i) {
            // Skip self (index 0 is the point itself when sorted)
            for (int j = 1; j < static_cast<int>(knn.indices[i].size()); ++j) {
                int nb = knn.indices[i][j];
                int lo = std::min(i, nb);
                int hi = std::max(i, nb);
                edge_set.insert({lo, hi});
            }
        }
    }

    // Build result with distances computed on normalized points
    AdjacencyResult result;
    result.edge_src.reserve(edge_set.size());
    result.edge_dst.reserve(edge_set.size());
    result.edge_dists.reserve(edge_set.size());

    for (auto& [lo, hi] : edge_set) {
        double dist = (norm_pts.col(lo) - norm_pts.col(hi)).norm();
        result.edge_src.push_back(lo);
        result.edge_dst.push_back(hi);
        result.edge_dists.push_back(dist);
    }

    if (filter) {
        filter_long_edges(result, n_mads);
    }

    return result;
}

} // namespace baysor
