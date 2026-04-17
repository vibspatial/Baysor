#pragma once

#include <string>
#include <stdexcept>

namespace baysor {

struct DataOptions {
    std::string x_col = "x";
    std::string y_col = "y";
    std::string z_col = "z";
    std::string gene_col = "gene";
    bool force_2d = false;
    int min_molecules_per_gene = 1;
    std::string exclude_genes;       // comma-separated, supports regex (e.g. "Blank*,MALAT1")
    int min_molecules_per_cell = 0;
    int min_molecules_per_segment = 0;
    int confidence_nn_id = 0;
};

struct SegmentationOptions {
    double scale = -1.0;
    std::string scale_std = "25%";
    bool estimate_scale_from_centers = true;
    int n_clusters = 4;
    double prior_segmentation_confidence = 0.2;
    int iters = 500;
    double tol = 0.0;     ///< convergence tolerance for main BMM loop (0 = run all iters, matching Julia)
    int n_cells_init = 0;
    std::string nuclei_genes;
    std::string cyto_genes;
    std::string unassigned_prior_label = "0";
};

struct PlottingOptions {
    int gene_composition_neighborhood = 0;
    int min_pixels_per_cell = 15;
    int max_plot_size = 3000;
    std::string ncv_method = "ri";
};

struct RunOptions {
    DataOptions data;
    SegmentationOptions segmentation;
    PlottingOptions plotting;
};

/// Parse scale_std: "25%" relative to scale, or absolute number like "5.0"
double parse_scale_std(const std::string& scale_std, double scale);

/// Fill in default parameter values derived from min_molecules_per_cell.
/// Throws if required info not available.
int default_param_value(const std::string& param, int min_molecules_per_cell,
                        int n_molecules = -1, int n_genes = -1);

/// Fill in zero-valued derived options from min_molecules_per_cell
void fill_and_check_data_options(DataOptions& opts);

/// Fill plotting defaults
void fill_and_check_plotting_options(PlottingOptions& opts, int min_molecules_per_cell,
                                     int n_genes = -1);

/// Load RunOptions from a TOML file. Missing keys retain defaults.
RunOptions load_config(const std::string& path);

/// Serialize RunOptions to a TOML dump file (like Julia's segmentation_params.dump.toml).
/// cli_cmd is the full CLI invocation string recorded as a comment.
void save_params_toml(const RunOptions& opts, const std::string& cli_cmd,
                      const std::string& path);

} // namespace baysor
