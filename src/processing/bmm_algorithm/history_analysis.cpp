#include "baysor/processing/bmm_algorithm/history_analysis.h"

namespace baysor {

template<int N>
std::pair<std::vector<int>, std::vector<double>> reassign_molecules_with_history(
    const BmmData<N>& data,
    int n_stable_iters,
    int min_molecules_per_cell,
    int max_iters,
    double outlier_confidence_threshold
) {
    // TODO: port from Julia history_analysis.jl
    return {{}, {}};
}

template std::pair<std::vector<int>, std::vector<double>>
    reassign_molecules_with_history<2>(const BmmData<2>&, int, int, int, double);
template std::pair<std::vector<int>, std::vector<double>>
    reassign_molecules_with_history<3>(const BmmData<3>&, int, int, int, double);

} // namespace baysor
