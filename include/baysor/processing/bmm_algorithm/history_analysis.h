#pragma once

#include "baysor/processing/models/bmm_data.h"
#include <vector>

namespace baysor {

/// Reassign molecules using frequency analysis of assignment history
template<int N>
std::pair<std::vector<int>, std::vector<double>> reassign_molecules_with_history(
    const BmmData<N>& data,
    int n_stable_iters,
    int min_molecules_per_cell = 10,
    int max_iters = 100,
    double outlier_confidence_threshold = 0.25
);

} // namespace baysor
