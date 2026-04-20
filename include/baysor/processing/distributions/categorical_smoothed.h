#pragma once

#include <vector>
#include <cmath>

namespace baysor {

/// Smoothed categorical distribution for gene expression profiles.
/// Each component (cell) maintains a count vector over genes.
struct CategoricalSmoothed {
    std::vector<int> gene_ids;   // sorted gene ids with explicit counts
    std::vector<float> counts;   // per-gene counts aligned with gene_ids
    float base_count = 0.0f;     // shared count for all genes (used for uniform init prior)
    double smooth = 1.0;
    double sum_counts = 0.0;
    int n_genes = 0;             // number of genes with nonzero counts
    int total_genes = 0;         // total number of genes in the panel

    explicit CategoricalSmoothed(int n_total_genes, double smooth = 1.0);

    double pdf(int gene_id, bool use_smoothing = true) const;
    float count_for_gene(int gene_id) const;
    int size() const { return total_genes; }

    void set_uniform_counts(float value);
    void set_dense_counts(const std::vector<float>& dense_counts);
    std::vector<float> dense_counts() const;

    /// M-step: recount from assigned molecule gene IDs
    void maximize(const int* gene_ids, int n, const double* confidences = nullptr);
    void maximize_indexed(const std::vector<int>& gene_ids, const int* ids, int n,
                          const std::vector<double>* confidences = nullptr);

    void reset();
};

} // namespace baysor
