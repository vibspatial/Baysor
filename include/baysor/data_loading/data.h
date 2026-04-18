#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <Eigen/Dense>
#include "baysor/utils/options.h"

namespace baysor {

/// Core molecule-level data. Struct-of-arrays layout for cache efficiency.
/// Replaces Julia's df_spatial DataFrame.
struct MoleculeData {
    // --- Spatial coordinates ---
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;   // empty if 2D

    // --- Gene identity ---
    std::vector<int> gene;                      // encoded as integer IDs (1-based, matching Julia convention)
    std::vector<std::string> gene_names;        // gene_names[gene[i]-1] gives the string name

    // --- Metadata columns (populated during processing) ---
    std::vector<double> confidence;
    std::vector<int> cluster;                    // MRF cluster assignment
    std::vector<int> prior_segmentation;         // prior segmentation labels (0 = unassigned)
    std::vector<double> nuclei_probs;

    // --- Extra string columns from input (preserved for output) ---
    // Column name -> values. Only populated if the input had extra columns we want to keep.
    // Omitted for now; can add later if needed.

    int n_molecules() const { return static_cast<int>(x.size()); }
    bool is_3d() const { return !z.empty(); }
    int n_dims() const { return is_3d() ? 3 : 2; }
    int n_genes() const { return static_cast<int>(gene_names.size()); }

    /// Return position data as Eigen matrix (dims x n_molecules), column-major.
    /// Each column is one molecule's coordinates.
    Eigen::MatrixXd position_matrix() const;
};

/// Load molecule data from CSV or Parquet file.
/// Handles column remapping, gene encoding, gene filtering, z-column handling.
/// This is the main entry point for data loading.
MoleculeData load_molecules(
    const std::string& path,
    const DataOptions& opts,
    const std::string& prior_column_name = "",
    const std::string& unassigned_prior_label = "0",
    int min_molecules_per_segment = 0
);

/// Internal: read a tabular file (CSV or Parquet) into raw column vectors.
/// Returns columns as named vectors. Detects format by file extension.
struct RawTableData {
    std::vector<double> x, y, z;       // z may be empty
    std::vector<std::string> gene_str; // raw gene strings
    // Prior segmentation column (if requested via :col_name syntax, loaded separately)
    bool has_z = false;
};
RawTableData read_tabular_file(const std::string& path, const DataOptions& opts);

/// Encode gene name strings to integer IDs (1-based).
/// Sorts unique names alphabetically, assigns 1..N.
/// Populates data.gene and data.gene_names.
void encode_genes(MoleculeData& data, const std::vector<std::string>& gene_strings);

/// Filter genes with fewer than min_molecules_per_gene total molecules.
/// Removes molecules of rare genes and re-encodes gene IDs.
void filter_genes_by_count(MoleculeData& data, int min_molecules_per_gene);

/// Filter genes matching any of the given patterns (supports * and ? wildcards).
/// Removes matching molecules and re-encodes gene IDs.
void filter_genes_by_pattern(MoleculeData& data, const std::vector<std::string>& patterns);

/// Read a string column from a CSV or Parquet file by name.
/// Used for loading prior segmentation from a column in the molecule file.
std::vector<std::string> read_string_column(const std::string& path, const std::string& col_name);

/// Read a numeric column from a CSV or Parquet file by name.
std::vector<double> read_double_column(const std::string& path, const std::string& col_name);

} // namespace baysor
