#include "baysor/processing/models/component.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace baysor {

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
        // Composition uses nuclei_probs as confidence weights (same as Julia)
        composition_params.maximize(gene_ids, n_points, nuclei_probs);
    }

    if (!freeze_position) {
        const ShapePrior<N>* prior_ptr = shape_prior.has_value() ? &shape_prior.value() : nullptr;
        position_params.maximize(pos_data, n_points, stride,
                                 nuclei_probs, prior_ptr, n_samples);
    }

    // Update component confidence from nuclei probabilities
    if (nuclei_probs != nullptr && n_points > 1) {
        std::vector<double> probs(nuclei_probs, nuclei_probs + n_points);
        // quantile at (1 - min_nuclei_frac)
        int k = static_cast<int>(std::floor((1.0 - min_nuclei_frac) * (n_points - 1)));
        k = std::max(0, std::min(k, n_points - 1));
        std::nth_element(probs.begin(), probs.begin() + k, probs.end());
        confidence = probs[k];
    }
}

template struct Component<2>;
template struct Component<3>;

} // namespace baysor
