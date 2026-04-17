#include "baysor/processing/models/component.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace baysor {

static double quantile_vec(std::vector<double> v, double p) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    double idx = p * (v.size() - 1);
    int lo = static_cast<int>(std::floor(idx));
    int hi = static_cast<int>(std::ceil(idx));
    if (lo == hi) return v[lo];
    double frac = idx - lo;
    return v[lo] * (1.0 - frac) + v[hi] * frac;
}

template<int N>
double Component<N>::pdf(const double* x, int gene, bool use_smoothing) const {
    double p = prior_probability * confidence * position_params.pdf(x);
    if (gene >= 0) {
        p *= composition_params.pdf(gene, use_smoothing);
    }
    return p;
}

template<int N>
double Component<N>::pdf_position_only(const double* x) const {
    return prior_probability * confidence * position_params.pdf(x);
}

template<int N>
void Component<N>::maximize(const double* pos_data, int stride,
                             const int* gene_ids, int n_points,
                             const double* nuclei_probs,
                             double min_nuclei_frac,
                             bool freeze_position,
                             bool freeze_composition) {
    n_samples = n_points;

    if (!freeze_composition) {
        // Julia maximizes composition with uniform weights (no nuclei_probs).
        // nuclei_probs are only used for position update and confidence below.
        composition_params.maximize(gene_ids, n_points, nullptr);
    }

    if (!freeze_position) {
        const ShapePrior<N>* prior_ptr = shape_prior.has_value() ? &shape_prior.value() : nullptr;
        position_params.maximize(pos_data, n_points, stride,
                                 nuclei_probs, prior_ptr, n_samples);
    }

    // Update component confidence from nuclei probabilities
    if (nuclei_probs != nullptr && n_points > 1) {
        std::vector<double> probs(nuclei_probs, nuclei_probs + n_points);
        confidence = quantile_vec(std::move(probs), 1.0 - min_nuclei_frac);
    }
}

template<int N>
void Component<N>::maximize_indexed(const Eigen::MatrixXd& pos_data,
                                    const std::vector<int>& gene_ids,
                                    const int* mol_ids, int n_points,
                                    const std::vector<double>* nuclei_probs,
                                    double min_nuclei_frac,
                                    bool freeze_position,
                                    bool freeze_composition) {
    n_samples = n_points;

    if (!freeze_composition) {
        composition_params.maximize_indexed(gene_ids, mol_ids, n_points, nullptr);
    }

    if (!freeze_position) {
        const ShapePrior<N>* prior_ptr = shape_prior.has_value() ? &shape_prior.value() : nullptr;
        position_params.maximize_indexed(pos_data, mol_ids, n_points,
                                         nuclei_probs, prior_ptr, n_samples);
    }

    if (nuclei_probs != nullptr && n_points > 1) {
        std::vector<double> probs;
        probs.reserve(n_points);
        for (int i = 0; i < n_points; ++i) {
            probs.push_back((*nuclei_probs)[mol_ids[i]]);
        }
        confidence = quantile_vec(std::move(probs), 1.0 - min_nuclei_frac);
    }
}

template struct Component<2>;
template struct Component<3>;

} // namespace baysor
