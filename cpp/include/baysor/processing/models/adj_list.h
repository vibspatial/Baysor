#pragma once

#include <vector>
#include <cstdint>

namespace baysor {

/// CSR-format adjacency list for the molecule graph (MRF).
///
/// Replaces Julia's Vector{Vector{Int}} with a flat, cache-friendly layout.
/// For molecule i, neighbors are indices[indptr[i]..indptr[i+1]),
/// with corresponding weights[indptr[i]..indptr[i+1]).
///
/// This layout is directly usable as a sparse matrix for GPU SpMM operations.
struct AdjList {
    std::vector<int32_t> indptr;    // size n_molecules + 1
    std::vector<int32_t> indices;   // size nnz (neighbor molecule IDs)
    std::vector<double> weights;    // size nnz (edge weights)

    int n_molecules() const { return static_cast<int>(indptr.size()) - 1; }
    int nnz() const { return static_cast<int>(indices.size()); }

    // Access neighbors of molecule i
    int neighbor_count(int i) const { return indptr[i + 1] - indptr[i]; }
    const int32_t* neighbor_ids(int i) const { return indices.data() + indptr[i]; }
    const double* neighbor_weights(int i) const { return weights.data() + indptr[i]; }

    /// Build from edge list (2 x n_edges) and edge weights
    static AdjList from_edge_list(
        const int* edges_src, const int* edges_dst,
        const double* edge_weights, int n_edges, int n_verts
    );
};

} // namespace baysor
