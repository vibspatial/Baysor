#include "baysor/processing/models/adj_list.h"
#include <algorithm>
#include <numeric>

namespace baysor {

AdjList AdjList::from_edge_list(
    const int* edges_src, const int* edges_dst,
    const double* edge_weights, int n_edges, int n_verts
) {
    if (n_edges == 0 || n_verts == 0) {
        AdjList adj;
        adj.indptr.assign(n_verts + 1, 0);
        return adj;
    }

    // Undirected: each edge (u, v) appears as both (u, v) and (v, u)
    int n_directed = n_edges * 2;

    struct DirEdge {
        int src, dst;
        double weight;
    };

    std::vector<DirEdge> edges(n_directed);
    for (int i = 0; i < n_edges; ++i) {
        edges[2 * i]     = {edges_src[i], edges_dst[i], edge_weights[i]};
        edges[2 * i + 1] = {edges_dst[i], edges_src[i], edge_weights[i]};
    }

    // Sort by (src, dst)
    std::sort(edges.begin(), edges.end(), [](const DirEdge& a, const DirEdge& b) {
        return a.src < b.src || (a.src == b.src && a.dst < b.dst);
    });

    // Deduplicate: same (src, dst) pair — keep the first occurrence
    auto last = std::unique(edges.begin(), edges.end(), [](const DirEdge& a, const DirEdge& b) {
        return a.src == b.src && a.dst == b.dst;
    });
    edges.erase(last, edges.end());
    n_directed = static_cast<int>(edges.size());

    // Build CSR
    AdjList adj;
    adj.indptr.assign(n_verts + 1, 0);
    adj.indices.resize(n_directed);
    adj.weights.resize(n_directed);

    // Count per row
    for (int i = 0; i < n_directed; ++i) {
        adj.indptr[edges[i].src + 1]++;
    }
    // Prefix sum
    for (int i = 1; i <= n_verts; ++i) {
        adj.indptr[i] += adj.indptr[i - 1];
    }
    // Fill
    for (int i = 0; i < n_directed; ++i) {
        adj.indices[i] = edges[i].dst;
        adj.weights[i] = edges[i].weight;
    }

    return adj;
}

} // namespace baysor
