#include "baysor/processing/distributions/categorical_smoothed.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace baysor {

namespace {

struct SparseCountScratch {
    std::vector<float> dense;
    std::vector<int> touched;

    void ensure_size(int n_total_genes) {
        if (static_cast<int>(dense.size()) < n_total_genes) {
            dense.resize(n_total_genes, 0.0f);
        }
        touched.clear();
    }
};

thread_local SparseCountScratch g_sparse_count_scratch;

template <typename GeneAccessor, typename WeightAccessor>
void rebuild_sparse_counts(
    CategoricalSmoothed& params,
    int n_points,
    GeneAccessor gene_at,
    WeightAccessor weight_at
) {
    auto& scratch = g_sparse_count_scratch;
    scratch.ensure_size(params.total_genes);

    params.gene_ids.clear();
    params.counts.clear();
    params.base_count = 0.0f;
    params.sum_counts = 0.0;
    params.n_genes = 0;

    for (int i = 0; i < n_points; ++i) {
        const int g = gene_at(i);
        if (g < 0 || g >= params.total_genes) continue;

        const float w = static_cast<float>(weight_at(i));
        if (scratch.dense[g] < 1e-10f) scratch.touched.push_back(g);
        scratch.dense[g] += w;
        params.sum_counts += static_cast<double>(w);
    }

    std::sort(scratch.touched.begin(), scratch.touched.end());
    params.gene_ids.reserve(scratch.touched.size());
    params.counts.reserve(scratch.touched.size());

    for (int g : scratch.touched) {
        const float cnt = scratch.dense[g];
        if (cnt >= 1e-10f) {
            params.gene_ids.push_back(g);
            params.counts.push_back(cnt);
        }
        scratch.dense[g] = 0.0f;
    }
    params.n_genes = static_cast<int>(params.gene_ids.size());
}

} // namespace

CategoricalSmoothed::CategoricalSmoothed(int n_total_genes, double smooth)
    : smooth(smooth), sum_counts(0.0), n_genes(0), total_genes(n_total_genes) {}

float CategoricalSmoothed::count_for_gene(int gene_id) const {
    if (gene_id < 0 || gene_id >= total_genes) return 0.0f;

    float cnt = base_count;
    auto it = std::lower_bound(gene_ids.begin(), gene_ids.end(), gene_id);
    if (it != gene_ids.end() && *it == gene_id) {
        cnt += counts[static_cast<size_t>(it - gene_ids.begin())];
    }
    return cnt;
}

void CategoricalSmoothed::set_uniform_counts(float value) {
    gene_ids.clear();
    counts.clear();
    base_count = value;
    if (value > 0.0f) {
        sum_counts = static_cast<double>(total_genes) * static_cast<double>(value);
        n_genes = total_genes;
    } else {
        sum_counts = 0.0;
        n_genes = 0;
    }
}

void CategoricalSmoothed::set_dense_counts(const std::vector<float>& dense_counts) {
    if (static_cast<int>(dense_counts.size()) != total_genes) {
        throw std::runtime_error("CategoricalSmoothed::set_dense_counts size mismatch");
    }

    gene_ids.clear();
    counts.clear();
    base_count = 0.0f;
    sum_counts = 0.0;
    n_genes = 0;

    for (int g = 0; g < total_genes; ++g) {
        const float cnt = dense_counts[g];
        if (cnt >= 1e-10f) {
            gene_ids.push_back(g);
            counts.push_back(cnt);
            ++n_genes;
            sum_counts += static_cast<double>(cnt);
        }
    }
}

std::vector<float> CategoricalSmoothed::dense_counts() const {
    std::vector<float> dense(static_cast<size_t>(total_genes), base_count);
    for (size_t i = 0; i < gene_ids.size(); ++i) {
        dense[static_cast<size_t>(gene_ids[i])] += counts[i];
    }
    return dense;
}

double CategoricalSmoothed::pdf(int gene_id, bool use_smoothing) const {
    if (gene_id < 0) return 1.0;  // missing gene — no composition term

    if (sum_counts < 1e-10) {
        // Uninformed prior: uniform over all genes
        return total_genes <= 0 ? 1.0 : 1.0 / static_cast<double>(total_genes);
    }

    double cnt = static_cast<double>(count_for_gene(gene_id));

    if (!use_smoothing) {
        return cnt / sum_counts;
    }

    // Laplace smoothing: only smooth up from 0 to `smooth`, don't add to already-nonzero
    // Matches Julia: max(cnt, smooth) / (sum_counts + smooth)
    return std::max(cnt, smooth) / (sum_counts + smooth);
}

void CategoricalSmoothed::maximize(const int* gene_ids, int n, const double* confidences) {
    rebuild_sparse_counts(
        *this, n,
        [&](int i) { return gene_ids[i]; },
        [&](int i) { return (confidences != nullptr) ? confidences[i] : 1.0; }
    );
}

void CategoricalSmoothed::maximize_indexed(const std::vector<int>& gene_ids, const int* ids, int n,
                                           const std::vector<double>* confidences) {
    rebuild_sparse_counts(
        *this, n,
        [&](int i) { return gene_ids[ids[i]]; },
        [&](int i) { return (confidences != nullptr) ? (*confidences)[ids[i]] : 1.0; }
    );
}

void CategoricalSmoothed::reset() {
    gene_ids.clear();
    counts.clear();
    base_count = 0.0f;
    sum_counts = 0.0;
    n_genes = 0;
}

} // namespace baysor
