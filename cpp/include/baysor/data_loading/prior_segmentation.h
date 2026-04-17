#pragma once

#include <string>
#include <vector>
#include <Eigen/Sparse>
#include "baysor/data_loading/data.h"

namespace baysor {

/// Determines the type of prior segmentation input
enum class PriorSegType {
    None,       // no prior segmentation
    Column,     // column in the molecule file (prefixed with ':')
    Image       // image file path (TIFF/PNG/MAT)
};

/// Detect prior segmentation type from the input string
PriorSegType detect_prior_seg_type(const std::string& prior_seg_arg);

/// Parse prior segmentation from a column in the molecule file.
/// The column values are strings or ints; unassigned_label marks unassigned molecules.
/// Returns 0-based integer labels (0 = unassigned).
/// Filters segments with fewer than min_molecules_per_segment molecules.
std::vector<int> parse_prior_from_column(
    const std::string& molecule_path,
    const std::string& col_name,
    const std::string& unassigned_label = "0",
    int min_molecules_per_segment = 0
);

/// Load prior segmentation from an image file (TIFF/PNG) or MAT file.
/// Returns (segment_per_molecule, optional_label_matrix).
/// Maps molecules to pixels via coordinate rounding.
/// Filters segments with fewer than min_molecules_per_segment molecules.
struct ImageSegResult {
    std::vector<int> segment_per_molecule;
    /// Pixel area (foreground pixel count) for each connected component,
    /// indexed by label-1 (index 0 -> label 1).
    /// Components filtered out by min_molecules_per_segment are zeroed.
    /// Empty for multi-label masks (not needed there).
    std::vector<size_t> component_pixel_areas;
};

/// Estimate scale from pixel areas of connected components (image-based prior).
/// Uses sqrt(area / π) as the radius of each nucleus.
/// Returns (median_radius, mad_radius).
std::pair<double, double> estimate_scale_from_image_areas(
    const std::vector<size_t>& component_pixel_areas
);

ImageSegResult load_prior_from_image(
    const std::string& image_path,
    const std::vector<double>& mol_x,
    const std::vector<double>& mol_y,
    int min_molecules_per_segment = 0
);

/// Top-level: load prior segmentation into MoleculeData.
/// Handles column-based and image-based priors.
/// Returns (scale, scale_std) if estimated from the prior, or (-1, -1) if not.
std::pair<double, double> load_prior_segmentation(
    MoleculeData& data,
    const std::string& prior_seg_arg,
    const std::string& molecule_path,
    const std::string& unassigned_label = "0",
    int min_molecules_per_segment = 0,
    int min_molecules_per_cell = 3,
    bool estimate_scale = true
);

/// Filter segments with too few molecules. Zeros out labels for small segments.
void filter_segmentation_labels(std::vector<int>& labels, int min_molecules_per_segment);

/// Estimate scale from inter-center distances of prior segments.
/// Returns (median_radius, mad_radius).
std::pair<double, double> estimate_scale_from_assignment(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& assignment,
    int min_molecules_per_cell
);

} // namespace baysor
