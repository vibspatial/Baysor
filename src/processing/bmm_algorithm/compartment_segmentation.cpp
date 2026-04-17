#include "baysor/processing/bmm_algorithm/compartment_segmentation.h"

namespace baysor {

std::pair<Eigen::MatrixXd, std::vector<bool>> init_nuclei_cyto_compartments(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const std::vector<std::string>& gene_names,
    const std::vector<std::string>& nuclei_genes,
    const std::vector<std::string>& cyto_genes,
    double scale
) {
    // TODO: port from Julia compartment_segmentation.jl
    return {};
}

CompartmentResult segment_molecule_compartments(
    Eigen::MatrixXd& assignment_probs,
    const std::vector<bool>& is_locked,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    double weight_mult,
    double tol,
    int max_iter,
    bool verbose
) {
    // TODO: port from Julia compartment_segmentation.jl
    return {};
}

void adjust_mrf_with_compartments(
    AdjList& adj_list,
    const std::vector<double>& nuclei_probs,
    const std::vector<double>& cyto_probs,
    double min_nuc_prob,
    double min_weight
) {
    // TODO: port from Julia compartment_segmentation.jl
}

} // namespace baysor
