#pragma once

#include "baysor/processing/models/adj_list.h"
#include <vector>
#include <Eigen/Dense>

namespace baysor {

/// Result of compartment (nuclei/cyto) segmentation
struct CompartmentResult {
    std::vector<int> assignment;         // 1=nuclei, 2=cyto, 3=unknown per molecule
    Eigen::MatrixXd assignment_probs;    // (n_compartments+1) x n_molecules
    std::vector<double> diffs;
    std::vector<double> change_fracs;
};

/// Initialize compartment assignment probabilities from nuclei/cyto gene markers
std::pair<Eigen::MatrixXd, std::vector<bool>> init_nuclei_cyto_compartments(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const std::vector<std::string>& gene_names,
    const std::vector<std::string>& nuclei_genes,
    const std::vector<std::string>& cyto_genes,
    double scale
);

/// Run MRF-based compartment segmentation
CompartmentResult segment_molecule_compartments(
    Eigen::MatrixXd& assignment_probs,
    const std::vector<bool>& is_locked,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    double weight_mult = 1.0,
    double tol = 0.01,
    int max_iter = 500,
    bool verbose = true
);

/// Adjust MRF edge weights based on nuclei/cyto probabilities
void adjust_mrf_with_compartments(
    AdjList& adj_list,
    const std::vector<double>& nuclei_probs,
    const std::vector<double>& cyto_probs,
    double min_nuc_prob = 0.25,
    double min_weight = 0.01
);

} // namespace baysor
