#include "baysor/processing/data_processing/boundary_estimation.h"

namespace baysor {

PolygonCollection boundary_polygons(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& cell_labels,
    const std::vector<std::string>* cell_names,
    double offset_rel
) {
    // TODO: port from Julia boundary_estimation.jl
    return {};
}

std::vector<Eigen::MatrixXd> boundary_polygons_from_grid(
    const Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic>& grid_labels,
    int grid_step
) {
    // TODO: port from Julia boundary_estimation.jl
    return {};
}

std::pair<PolygonCollection, PolygonCollection> boundary_polygons_auto(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& assignment,
    bool estimate_per_z,
    const std::vector<std::string>* cell_names,
    bool verbose
) {
    // TODO: port from Julia boundary_estimation.jl
    return {};
}

} // namespace baysor
