#pragma once

#include "baysor/processing/data_processing/boundary_estimation.h"
#include <Eigen/Sparse>
#include <map>
#include <string>
#include <variant>
#include <vector>

namespace baysor {

struct MoleculeData;

/// Output file paths (mirrors Julia OutputPaths)
struct OutputPaths {
    std::string segmented_df;
    std::string cell_stats;
    std::string counts;
    std::string diagnostic_report;
    std::string molecule_plot;
    std::string polygons_2d;
    std::string polygons_3d;
    std::string params_dump;   ///< segmentation_params.dump.toml
    std::string log_file;      ///< segmentation_log.log
};

OutputPaths get_output_paths(const std::string& base_path, const std::string& count_matrix_format);

/// Save segmented molecule table to CSV.
/// Optional columns: ncv_color (hex strings), assignment_confidence (0..1),
/// cluster (1-based integer cluster ID).
void save_segmented_df(const MoleculeData& data,
                       const std::vector<int>& assignment,
                       const std::vector<std::string>& gene_names,
                       const std::string& path,
                       const std::vector<std::string>* ncv_color = nullptr,
                       const std::vector<double>* assignment_confidence = nullptr,
                       const std::vector<int>* cluster = nullptr);

/// Save cell statistics to CSV
void save_cell_stat_df(const Eigen::MatrixXd& stats,
                       const std::vector<std::string>& cell_names,
                       const std::vector<std::string>& col_names,
                       const std::string& path);

/// Per-column attributes that can be written into a Loom col_attrs group.
/// Each value is either a vector of strings or a vector of doubles.
using LoomColAttrs = std::map<std::string,
    std::variant<std::vector<std::string>, std::vector<double>>>;

/// Save count matrix to Loom (HDF5) format.
/// matrix  : n_cells × n_genes  (rows = cells/NCVs, cols = genes)
/// col_attrs: optional extra per-cell attributes (e.g. ncv_color, confidence)
void save_matrix_to_loom(const Eigen::SparseMatrix<float>& matrix,
                         const std::vector<std::string>& gene_names,
                         const std::vector<std::string>& cell_names,
                         const std::string& path,
                         const LoomColAttrs& col_attrs = {});

/// Save count matrix to TSV
void save_matrix_to_tsv(const Eigen::SparseMatrix<double>& matrix,
                        const std::vector<std::string>& gene_names,
                        const std::vector<std::string>& cell_names,
                        const std::string& path);

/// Convert polygons to GeoJSON and save
void save_polygons_geojson(const PolygonCollection& polygons,
                           const std::string& path,
                           const std::string& format = "FeatureCollection");

} // namespace baysor
