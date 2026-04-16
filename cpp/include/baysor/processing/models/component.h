#pragma once

#include "baysor/processing/distributions/mv_normal.h"
#include "baysor/processing/distributions/categorical_smoothed.h"
#include <optional>
#include <unordered_map>

namespace baysor {

/// A cell component: spatial distribution + gene expression profile.
/// Template parameter N = spatial dimensionality (2 or 3).
template<int N>
struct Component {
    MvNormal<N> position_params;
    CategoricalSmoothed composition_params;

    int n_samples = 0;
    double prior_probability = 1.0;
    double confidence = 1.0;
    int guid = -1;

    /// Owned shape prior for covariance regularization (optional).
    /// Stored by value to avoid pointer lifetime issues.
    std::optional<ShapePrior<N>> shape_prior = std::nullopt;

    // Per-segment molecule counts (for prior segmentation adjustment)
    std::unordered_map<int, int> n_molecules_per_segment;

    Component(const MvNormal<N>& pos, const CategoricalSmoothed& comp,
              std::optional<ShapePrior<N>> prior = std::nullopt, int guid = -1)
        : position_params(pos), composition_params(comp), shape_prior(std::move(prior)), guid(guid) {}

    /// Likelihood of observing molecule at position x with gene g
    double pdf(const double* x, int gene, bool use_smoothing = true) const;
    /// Likelihood with missing gene (position only)
    double pdf_position_only(const double* x) const;

    /// M-step: re-estimate from assigned molecules
    void maximize(const double* pos_data, int stride,
                  const int* gene_ids, int n_points,
                  const double* nuclei_probs = nullptr,
                  double min_nuclei_frac = 0.1,
                  bool freeze_position = false,
                  bool freeze_composition = false);
};

extern template struct Component<2>;
extern template struct Component<3>;

} // namespace baysor
