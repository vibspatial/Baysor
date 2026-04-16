#pragma once

#include <vector>
#include <cmath>

namespace baysor {

/// Smoothed categorical distribution for gene expression profiles.
/// Each component (cell) maintains a count vector over genes.
struct CategoricalSmoothed {
    std::vector<double> counts;  // per-gene counts (float to allow weighted)
    double smooth = 1.0;
    double sum_counts = 0.0;
    int n_genes = 0;             // number of genes with nonzero counts

    explicit CategoricalSmoothed(int n_total_genes, double smooth = 1.0);

    double pdf(int gene_id, bool use_smoothing = true) const;

    /// M-step: recount from assigned molecule gene IDs
    void maximize(const int* gene_ids, int n, const double* confidences = nullptr);

    void reset();
};

} // namespace baysor
