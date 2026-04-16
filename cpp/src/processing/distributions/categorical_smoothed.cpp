#include "baysor/processing/distributions/categorical_smoothed.h"
#include <algorithm>
#include <cmath>

namespace baysor {

CategoricalSmoothed::CategoricalSmoothed(int n_total_genes, double smooth)
    : counts(n_total_genes, 0.0), smooth(smooth), sum_counts(0.0), n_genes(0) {}

double CategoricalSmoothed::pdf(int gene_id, bool use_smoothing) const {
    if (gene_id < 0) return 1.0;  // missing gene — no composition term

    if (sum_counts < 1e-10) {
        // Uninformed prior: uniform over all genes
        return counts.empty() ? 1.0 : 1.0 / static_cast<double>(counts.size());
    }

    double cnt = (gene_id < static_cast<int>(counts.size())) ? counts[gene_id] : 0.0;

    if (!use_smoothing) {
        return cnt / sum_counts;
    }

    // Laplace smoothing: only smooth up from 0 to `smooth`, don't add to already-nonzero
    // Matches Julia: max(cnt, smooth) / (sum_counts + smooth)
    return std::max(cnt, smooth) / (sum_counts + smooth);
}

void CategoricalSmoothed::maximize(const int* gene_ids, int n, const double* confidences) {
    // Reset
    std::fill(counts.begin(), counts.end(), 0.0);
    sum_counts = 0.0;
    n_genes = 0;

    for (int i = 0; i < n; ++i) {
        int g = gene_ids[i];
        if (g < 0 || g >= static_cast<int>(counts.size())) continue;  // missing or out-of-range

        double w = (confidences != nullptr) ? confidences[i] : 1.0;
        if (counts[g] < 1e-10) n_genes++;  // first nonzero entry for this gene
        counts[g] += w;
        sum_counts += w;
    }
}

void CategoricalSmoothed::reset() {
    std::fill(counts.begin(), counts.end(), 0.0);
    sum_counts = 0.0;
    n_genes = 0;
}

} // namespace baysor
