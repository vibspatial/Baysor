#include <gtest/gtest.h>

#include "baysor/reporting/output.h"
#include "baysor/processing/data_processing/triangulation.h"
#include "baysor/segmentation/segmentation.h"
#include "baysor/utils/general.h"
#include "baysor/utils/options.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef BAYSOR_TEST_CLI_PATH
#error "BAYSOR_TEST_CLI_PATH must identify the Baysor CLI used by the baseline test"
#endif

#ifndef BAYSOR_TEST_FIXTURE_DIR
#error "BAYSOR_TEST_FIXTURE_DIR must identify the native baseline fixture"
#endif

#ifndef BAYSOR_TEST_SCIENTIFIC_PARITY_DIR
#error "BAYSOR_TEST_SCIENTIFIC_PARITY_DIR must identify the scientific parity fixtures"
#endif

#ifndef BAYSOR_TEST_CMAKE_COMMAND
#error "BAYSOR_TEST_CMAKE_COMMAND must identify CMake for portable environment setup"
#endif

namespace {

namespace fs = std::filesystem;

struct DelimitedTable {
    std::vector<std::string> columns;
    std::unordered_map<std::string, std::size_t> column_index;
    std::vector<std::vector<std::string>> rows;

    const std::string& value(const std::vector<std::string>& row, const std::string& column) const {
        auto it = column_index.find(column);
        if (it == column_index.end()) {
            throw std::runtime_error("Missing column '" + column + "'");
        }
        return row.at(it->second);
    }
};

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::vector<std::string> split_delimited(const std::string& line, char delimiter) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, delimiter)) fields.push_back(field);
    if (!line.empty() && line.back() == delimiter) fields.emplace_back();
    return fields;
}

DelimitedTable read_delimited_table(const fs::path& path, char delimiter) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open table: " + path.string());

    DelimitedTable table;
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("Table is empty: " + path.string());
    }
    if (!line.empty() && line.back() == '\r') line.pop_back();
    table.columns = split_delimited(line, delimiter);
    for (std::size_t i = 0; i < table.columns.size(); ++i) {
        require(table.column_index.emplace(table.columns[i], i).second,
                "Duplicate column '" + table.columns[i] + "' in " + path.string());
    }

    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        auto row = split_delimited(line, delimiter);
        require(row.size() == table.columns.size(),
                "Unexpected column count in " + path.string() + ": " + line);
        table.rows.push_back(std::move(row));
    }
    return table;
}

std::map<std::string, const std::vector<std::string>*> rows_by_key(
    const DelimitedTable& table,
    const std::string& key_column
) {
    std::map<std::string, const std::vector<std::string>*> result;
    for (const auto& row : table.rows) {
        const auto& key = table.value(row, key_column);
        require(result.emplace(key, &row).second, "Duplicate key '" + key + "'");
    }
    return result;
}

double parse_number(const std::string& value, const std::string& context) {
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    // strtod preserves a representable subnormal result even when the C library
    // reports range loss; those tiny confidence values are valid test data.
    require(end != value.c_str() && end == value.c_str() + value.size(),
            "Invalid number '" + value + "' for " + context);
    return parsed;
}

bool parse_boolean(const std::string& value, const std::string& context) {
    if (value == "true" || value == "1") return true;
    if (value == "false" || value == "0") return false;
    throw std::runtime_error("Invalid boolean '" + value + "' for " + context);
}

void require_near(
    double actual,
    double expected,
    double absolute_tolerance,
    double relative_tolerance,
    const std::string& context
) {
    const double tolerance = absolute_tolerance + relative_tolerance * std::abs(expected);
    require(std::abs(actual - expected) <= tolerance,
            context + " differs: expected " + std::to_string(expected) +
                ", got " + std::to_string(actual));
}

using CellLabelMap = std::map<std::string, std::string>;  // candidate -> reference

void record_bijective_label(
    CellLabelMap& labels,
    CellLabelMap& reverse_labels,
    const std::string& actual_label,
    const std::string& reference_label,
    const std::string& context
) {
    auto [mapping, inserted] = labels.emplace(actual_label, reference_label);
    require(inserted || mapping->second == reference_label,
            "Candidate " + context + " " + actual_label + " combines reference partitions");
    auto [reverse, reverse_inserted] = reverse_labels.emplace(reference_label, actual_label);
    require(reverse_inserted || reverse->second == actual_label,
            "Reference " + context + " " + reference_label + " was split across candidates");
}

CellLabelMap compare_molecule_tables(
    const fs::path& actual_path,
    const fs::path& reference_path,
    CellLabelMap* cluster_labels = nullptr
) {
    const auto actual = read_delimited_table(actual_path, ',');
    const auto reference = read_delimited_table(reference_path, ',');
    require(actual.columns == reference.columns, "Segmented-molecule columns changed");

    const auto actual_rows = rows_by_key(actual, "transcript_id");
    const auto reference_rows = rows_by_key(reference, "transcript_id");
    require(actual_rows.size() == reference_rows.size(), "Retained transcript count changed");

    CellLabelMap labels;
    std::map<std::string, std::string> reverse_labels;
    CellLabelMap clusters;
    CellLabelMap reverse_clusters;
    for (const auto& [transcript_id, reference_row] : reference_rows) {
        auto actual_it = actual_rows.find(transcript_id);
        require(actual_it != actual_rows.end(), "Missing transcript_id " + transcript_id);
        const auto& actual_row = *actual_it->second;

        require(actual.value(actual_row, "gene") == reference.value(*reference_row, "gene"),
                "Gene changed for transcript_id " + transcript_id);
        std::vector<std::string> coordinates = {"x", "y"};
        if (actual.column_index.count("z") != 0) coordinates.push_back("z");
        for (const auto& coordinate : coordinates) {
            require_near(
                parse_number(actual.value(actual_row, coordinate), coordinate),
                parse_number(reference.value(*reference_row, coordinate), coordinate),
                1e-9,
                1e-9,
                std::string(coordinate) + " for transcript_id " + transcript_id
            );
        }
        if (actual.column_index.count("cluster") != 0) {
            record_bijective_label(
                clusters,
                reverse_clusters,
                actual.value(actual_row, "cluster"),
                reference.value(*reference_row, "cluster"),
                "molecule cluster"
            );
        }
        for (const auto& confidence : {"confidence", "assignment_confidence"}) {
            require_near(
                parse_number(actual.value(actual_row, confidence), confidence),
                parse_number(reference.value(*reference_row, confidence), confidence),
                1e-8,
                1e-7,
                std::string(confidence) + " for transcript_id " + transcript_id
            );
        }

        const auto& actual_cell = actual.value(actual_row, "cell");
        const auto& reference_cell = reference.value(*reference_row, "cell");
        const bool actual_noise = parse_boolean(actual.value(actual_row, "is_noise"), "is_noise");
        const bool reference_noise = parse_boolean(reference.value(*reference_row, "is_noise"), "is_noise");
        require(actual_noise == reference_noise, "Noise assignment changed for transcript_id " + transcript_id);
        require(actual_noise == (actual_cell == "0"),
                "Candidate cell and is_noise disagree for transcript_id " + transcript_id);
        require(reference_noise == (reference_cell == "0"),
                "Reference cell and is_noise disagree for transcript_id " + transcript_id);

        if (!reference_noise) {
            record_bijective_label(
                labels, reverse_labels, actual_cell, reference_cell, "cell");
        }
    }
    require(labels.size() == reverse_labels.size(), "Cell relabelling is not bijective");
    require(clusters.size() == reverse_clusters.size(), "Cluster relabelling is not bijective");
    if (cluster_labels != nullptr) *cluster_labels = std::move(clusters);
    return labels;
}

void compare_cell_statistics(
    const fs::path& actual_path,
    const fs::path& reference_path,
    const CellLabelMap& labels,
    const CellLabelMap* cluster_labels = nullptr
) {
    const auto actual = read_delimited_table(actual_path, ',');
    const auto reference = read_delimited_table(reference_path, ',');
    require(actual.columns == reference.columns, "Cell-statistics columns changed");

    const auto actual_rows = rows_by_key(actual, "cell");
    const auto reference_rows = rows_by_key(reference, "cell");
    require(actual_rows.size() == reference_rows.size(), "Cell-statistics row count changed");

    for (const auto& [actual_cell, actual_row] : actual_rows) {
        auto label_it = labels.find(actual_cell);
        require(label_it != labels.end(), "Cell statistics contain unknown cell " + actual_cell);
        auto reference_it = reference_rows.find(label_it->second);
        require(reference_it != reference_rows.end(), "Reference cell statistics are incomplete");

        for (const auto& column : actual.columns) {
            if (column == "cell") continue;
            if (column == "cluster") {
                require(cluster_labels != nullptr,
                        "Cell statistics contain clusters without a molecule-cluster mapping");
                const auto actual_cluster = actual.value(*actual_row, column);
                auto cluster_it = cluster_labels->find(actual_cluster);
                require(cluster_it != cluster_labels->end(),
                        "Cell statistics contain unknown cluster " + actual_cluster);
                require(cluster_it->second == reference.value(*reference_it->second, column),
                        "Cell cluster changed for cell " + actual_cell);
                continue;
            }
            require_near(
                parse_number(actual.value(*actual_row, column), column),
                parse_number(reference.value(*reference_it->second, column), column),
                1e-5,
                1e-6,
                std::string(column) + " for cell " + actual_cell
            );
        }
    }
}

void compare_count_matrices(
    const fs::path& actual_path,
    const fs::path& reference_path,
    const CellLabelMap& labels
) {
    const auto actual = read_delimited_table(actual_path, '\t');
    const auto reference = read_delimited_table(reference_path, '\t');
    require(actual.columns.size() == reference.columns.size(), "Count-matrix cell count changed");
    require(actual.columns.front() == "gene" && reference.columns.front() == "gene",
            "Count matrix does not start with the gene column");

    const auto actual_rows = rows_by_key(actual, "gene");
    const auto reference_rows = rows_by_key(reference, "gene");
    require(actual_rows.size() == reference_rows.size(), "Count-matrix gene count changed");

    for (const auto& [gene, reference_row] : reference_rows) {
        auto actual_it = actual_rows.find(gene);
        require(actual_it != actual_rows.end(), "Count matrix is missing gene " + gene);
        for (std::size_t column = 1; column < actual.columns.size(); ++column) {
            const auto& actual_cell = actual.columns[column];
            auto label_it = labels.find(actual_cell);
            require(label_it != labels.end(), "Count matrix contains unknown cell " + actual_cell);
            auto reference_column = reference.column_index.find(label_it->second);
            require(reference_column != reference.column_index.end(),
                    "Reference count matrix is missing cell " + label_it->second);
            require_near(
                parse_number(actual_it->second->at(column), "count"),
                parse_number(reference_row->at(reference_column->second), "count"),
                0.0,
                0.0,
                "count for gene " + gene + " and cell " + actual_cell
            );
        }
    }
}

using Point = std::array<double, 2>;
using PolygonMap = std::map<std::string, std::vector<Point>>;

PolygonMap read_polygons(const nlohmann::json& document) {
    require(document.at("type") == "FeatureCollection", "Expected a GeoJSON FeatureCollection");

    PolygonMap polygons;
    for (const auto& feature : document.at("features")) {
        const std::string cell = feature.at("properties").at("cell").get<std::string>();
        const auto& ring = feature.at("geometry").at("coordinates").at(0);
        std::vector<Point> vertices;
        vertices.reserve(ring.size());
        for (const auto& coordinate : ring) {
            vertices.push_back({coordinate.at(0).get<double>(), coordinate.at(1).get<double>()});
        }
        if (vertices.size() > 1 && vertices.front() == vertices.back()) vertices.pop_back();
        std::sort(vertices.begin(), vertices.end());
        require(polygons.emplace(cell, std::move(vertices)).second,
                "Duplicate polygon for cell " + cell);
    }
    return polygons;
}

PolygonMap read_polygons(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open polygons: " + path.string());
    nlohmann::json document;
    input >> document;
    return read_polygons(document);
}

void compare_polygon_maps(
    const PolygonMap& actual,
    const PolygonMap& reference,
    const CellLabelMap& labels,
    const std::string& context
) {
    require(actual.size() == reference.size(), context + " polygon count changed");

    for (const auto& [actual_cell, actual_vertices] : actual) {
        auto label_it = labels.find(actual_cell);
        require(label_it != labels.end(), context + " contain unknown cell " + actual_cell);
        auto reference_it = reference.find(label_it->second);
        require(reference_it != reference.end(), "Reference " + context + " are missing a cell");
        const auto& reference_vertices = reference_it->second;
        require(actual_vertices.size() == reference_vertices.size(),
                context + " vertex count changed for cell " + actual_cell);
        for (std::size_t i = 0; i < actual_vertices.size(); ++i) {
            for (std::size_t dimension = 0; dimension < 2; ++dimension) {
                require_near(
                    actual_vertices[i][dimension],
                    reference_vertices[i][dimension],
                    1e-8,
                    1e-8,
                    context + " coordinate for cell " + actual_cell
                );
            }
        }
    }
}

void compare_polygons(
    const fs::path& actual_path,
    const fs::path& reference_path,
    const CellLabelMap& labels
) {
    const auto actual = read_polygons(actual_path);
    const auto reference = read_polygons(reference_path);
    compare_polygon_maps(actual, reference, labels, "Boundaries");
}

void compare_polygon_stacks(
    const fs::path& actual_path,
    const fs::path& reference_path,
    const CellLabelMap& labels
) {
    std::ifstream actual_input(actual_path);
    std::ifstream reference_input(reference_path);
    if (!actual_input) throw std::runtime_error("Cannot open polygon stack: " + actual_path.string());
    if (!reference_input) {
        throw std::runtime_error("Cannot open polygon stack: " + reference_path.string());
    }
    nlohmann::json actual;
    nlohmann::json reference;
    actual_input >> actual;
    reference_input >> reference;
    require(actual.is_object() && reference.is_object(), "Expected object-valued polygon stacks");
    require(actual.size() == reference.size(), "3D boundary layer count changed");
    for (auto reference_layer = reference.begin(); reference_layer != reference.end(); ++reference_layer) {
        auto actual_layer = actual.find(reference_layer.key());
        require(actual_layer != actual.end(), "Missing 3D boundary layer " + reference_layer.key());
        compare_polygon_maps(
            read_polygons(*actual_layer),
            read_polygons(reference_layer.value()),
            labels,
            "3D boundaries in layer " + reference_layer.key()
        );
    }
}

void compare_resolved_options(const fs::path& actual_path, const fs::path& reference_path) {
    const baysor::RunOptions actual = baysor::load_config(actual_path.string());
    const baysor::RunOptions reference = baysor::load_config(reference_path.string());

    require(actual.molecules.min_molecules_per_gene == reference.molecules.min_molecules_per_gene,
            "Resolved min_molecules_per_gene changed");
    require(actual.molecules.min_molecules_per_cell == reference.molecules.min_molecules_per_cell,
            "Resolved min_molecules_per_cell changed");
    require(actual.molecules.confidence_nn_id == reference.molecules.confidence_nn_id,
            "Resolved confidence_nn_id changed");
    require(actual.prior.type == reference.prior.type, "Resolved prior type changed");
    require(actual.prior.column_name == reference.prior.column_name,
            "Resolved prior column changed");
    require(actual.prior.unassigned_label == reference.prior.unassigned_label,
            "Resolved unassigned prior label changed");
    require(actual.prior.min_molecules_per_segment == reference.prior.min_molecules_per_segment,
            "Resolved min_molecules_per_segment changed");
    require(actual.prior.estimate_scale_from_prior == reference.prior.estimate_scale_from_prior,
            "Resolved prior scale policy changed");
    require_near(actual.segmentation.scale, reference.segmentation.scale, 0.0, 0.0,
                 "Resolved scale");
    require(actual.segmentation.scale_std == reference.segmentation.scale_std,
            "Resolved scale_std changed");
    require(actual.segmentation.cluster_method == reference.segmentation.cluster_method,
            "Resolved cluster method changed");
    require(actual.segmentation.n_clusters == reference.segmentation.n_clusters,
            "Resolved cluster count changed");
    require_near(actual.segmentation.cluster_resolution,
                 reference.segmentation.cluster_resolution, 0.0, 0.0,
                 "Resolved cluster resolution");
    require(actual.segmentation.cluster_graph_k == reference.segmentation.cluster_graph_k,
            "Resolved cluster graph k changed");
    require(actual.segmentation.cluster_n_dims == reference.segmentation.cluster_n_dims,
            "Resolved cluster dimension count changed");
    require(actual.segmentation.cluster_basis_sample_size ==
                reference.segmentation.cluster_basis_sample_size,
            "Resolved cluster basis sample size changed");
    require_near(
        actual.segmentation.prior_segmentation_confidence,
        reference.segmentation.prior_segmentation_confidence,
        0.0,
        0.0,
        "Resolved prior segmentation confidence"
    );
    require(actual.segmentation.iters == reference.segmentation.iters,
            "Resolved iteration count changed");
    require(actual.segmentation.n_cells_init == reference.segmentation.n_cells_init,
            "Resolved n_cells_init changed");
    require(actual.plotting.gene_composition_neighborhood ==
                reference.plotting.gene_composition_neighborhood,
            "Resolved gene composition neighborhood changed");
}

void compare_resolved_options(
    const baysor::ResolvedSegmentationOptions& actual,
    const fs::path& reference_path
) {
    const baysor::RunOptions reference = baysor::load_config(reference_path.string());

    require(actual.molecules.min_molecules_per_gene == reference.molecules.min_molecules_per_gene,
            "Resolved min_molecules_per_gene changed");
    require(actual.molecules.min_molecules_per_cell == reference.molecules.min_molecules_per_cell,
            "Resolved min_molecules_per_cell changed");
    require(actual.molecules.confidence_nn_id == reference.molecules.confidence_nn_id,
            "Resolved confidence_nn_id changed");
    require(actual.prior.type == reference.prior.type, "Resolved prior type changed");
    require(actual.prior.column_name == reference.prior.column_name,
            "Resolved prior column changed");
    require(actual.prior.unassigned_label == reference.prior.unassigned_label,
            "Resolved unassigned prior label changed");
    require(actual.prior.min_molecules_per_segment == reference.prior.min_molecules_per_segment,
            "Resolved min_molecules_per_segment changed");
    require(actual.prior.estimate_scale_from_prior == reference.prior.estimate_scale_from_prior,
            "Resolved prior scale policy changed");
    require_near(actual.segmentation.scale, reference.segmentation.scale, 0.0, 0.0,
                 "Resolved scale");
    require(actual.segmentation.scale_std == reference.segmentation.scale_std,
            "Resolved scale_std changed");
    require(actual.segmentation.cluster_method == reference.segmentation.cluster_method,
            "Resolved cluster method changed");
    require(actual.segmentation.n_clusters == reference.segmentation.n_clusters,
            "Resolved cluster count changed");
    require_near(actual.segmentation.cluster_resolution,
                 reference.segmentation.cluster_resolution, 0.0, 0.0,
                 "Resolved cluster resolution");
    require(actual.segmentation.cluster_graph_k == reference.segmentation.cluster_graph_k,
            "Resolved cluster graph k changed");
    require(actual.segmentation.cluster_n_dims == reference.segmentation.cluster_n_dims,
            "Resolved cluster dimension count changed");
    require(actual.segmentation.cluster_basis_sample_size ==
                reference.segmentation.cluster_basis_sample_size,
            "Resolved cluster basis sample size changed");
    require_near(
        actual.segmentation.prior_segmentation_confidence,
        reference.segmentation.prior_segmentation_confidence,
        0.0,
        0.0,
        "Resolved prior segmentation confidence");
    require(actual.segmentation.iters == reference.segmentation.iters,
            "Resolved iteration count changed");
    require(actual.segmentation.n_cells_init == reference.segmentation.n_cells_init,
            "Resolved n_cells_init changed");
    require(actual.neighborhood_composition.neighborhood_size ==
                reference.plotting.gene_composition_neighborhood,
            "Resolved gene composition neighborhood changed");
}

std::string quote_argument(const fs::path& value) {
    const std::string text = value.string();
    require(text.find('"') == std::string::npos, "Test path contains an unsupported quote: " + text);
    return "\"" + text + "\"";
}

class TemporaryDirectory {
public:
    TemporaryDirectory() {
        const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = fs::temp_directory_path() / ("baysor_native_baseline_" + std::to_string(nonce));
        require(fs::create_directories(path_), "Cannot create test output directory: " + path_.string());
    }

    ~TemporaryDirectory() {
        std::error_code error;
        fs::remove_all(path_, error);
    }

    const fs::path& path() const { return path_; }

private:
    fs::path path_;
};

TEST(NativeBaseline, PreExtractionCliMatchesSemanticReference) {
    const fs::path fixture_dir(BAYSOR_TEST_FIXTURE_DIR);
    const fs::path reference_dir = fixture_dir / "reference";
    TemporaryDirectory output;

    const std::string command =
        quote_argument(fs::path(BAYSOR_TEST_CMAKE_COMMAND)) +
        " -E env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE " +
        quote_argument(fs::path(BAYSOR_TEST_CLI_PATH)) +
        " run " + quote_argument(fixture_dir / "molecules.csv") +
        " --config " + quote_argument(fixture_dir / "config.toml") +
        " --output " + quote_argument(output.path()) +
        " --output-style legacy --polygon-format FeatureCollection" +
        " --count-matrix-format tsv --skip-ncv-color";

    ASSERT_EQ(std::system(command.c_str()), 0) << "Baseline CLI command failed: " << command;

    try {
        const auto labels = compare_molecule_tables(
            output.path() / "segmentation.csv",
            reference_dir / "segmentation.csv"
        );
        compare_cell_statistics(
            output.path() / "segmentation_cell_stats.csv",
            reference_dir / "segmentation_cell_stats.csv",
            labels
        );
        compare_count_matrices(
            output.path() / "segmentation_counts.tsv",
            reference_dir / "segmentation_counts.tsv",
            labels
        );
        compare_polygons(
            output.path() / "segmentation_polygons_2d.json",
            reference_dir / "segmentation_polygons_2d.json",
            labels
        );
        compare_resolved_options(
            output.path() / "segmentation_params.dump.toml",
            reference_dir / "resolved_params.toml"
        );
    } catch (const std::exception& error) {
        FAIL() << error.what();
    }
}

TEST(SegmentationOperation, DirectCallMatchesNativeBaseline) {
    const fs::path fixture_dir(BAYSOR_TEST_FIXTURE_DIR);
    const fs::path reference_dir = fixture_dir / "reference";
    const baysor::RunOptions config = baysor::load_config(
        (fixture_dir / "config.toml").string());

    baysor::SegmentationRequest request;
    request.molecules.path = (fixture_dir / "molecules.csv").string();
    request.molecules.options = config.molecules;
    request.prior = config.prior;
    request.segmentation = config.segmentation;
    request.neighborhood_composition.neighborhood_size =
        config.plotting.gene_composition_neighborhood;
    request.neighborhood_composition.method = config.plotting.ncv_method;
    request.requested_products = baysor::SegmentationProducts::none();
    request.requested_products.molecule_assignments = true;
    request.requested_products.molecule_confidence = true;
    request.requested_products.assignment_confidence = true;
    request.requested_products.cell_statistics = true;
    request.requested_products.boundaries = true;
    request.requested_products.count_matrix = true;
    request.requested_products.diagnostics = true;
    request.random_seed = baysor::kDefaultSegmentationSeed;
    request.execution.native_threads = 1;

    const auto outcome = baysor::run_segmentation(request, baysor::CancellationToken{});
    ASSERT_TRUE(std::holds_alternative<baysor::SegmentationResult>(outcome));
    const auto& result = std::get<baysor::SegmentationResult>(outcome);

    EXPECT_EQ(result.molecules.n_molecules(), static_cast<int>(result.cell_assignments.size()));
    EXPECT_EQ(result.molecules.n_molecules(), static_cast<int>(result.molecules.confidence.size()));
    EXPECT_EQ(result.molecules.n_molecules(), static_cast<int>(result.assignment_confidence.size()));
    EXPECT_TRUE(result.cell_statistics.has_value());
    EXPECT_TRUE(result.boundaries_2d.has_value());
    EXPECT_FALSE(result.boundaries_3d.has_value());
    EXPECT_TRUE(result.count_matrix.has_value());
    EXPECT_TRUE(result.diagnostics.has_value());
    EXPECT_TRUE(result.produced_products.molecule_assignments);
    EXPECT_TRUE(result.produced_products.molecule_confidence);
    EXPECT_TRUE(result.produced_products.assignment_confidence);
    EXPECT_FALSE(result.produced_products.molecule_clusters);
    EXPECT_FALSE(result.produced_products.neighborhood_composition_colors);
    EXPECT_TRUE(result.produced_products.cell_statistics);
    EXPECT_TRUE(result.produced_products.boundaries);
    EXPECT_TRUE(result.produced_products.count_matrix);
    EXPECT_TRUE(result.produced_products.diagnostics);
    EXPECT_EQ(result.provenance.random_seed, baysor::kDefaultSegmentationSeed);
    EXPECT_EQ(result.provenance.random_substream_contract_version,
              baysor::kRandomSubstreamContractVersion);
    EXPECT_EQ(result.provenance.random_substreams.size(), 4U);
    EXPECT_EQ(result.provenance.effective_native_threads, 1);
    EXPECT_EQ(result.provenance.baysor_version, "0.8.3");

    TemporaryDirectory output;
    ASSERT_NO_THROW(baysor::save_segmented_df(
        result, (output.path() / "segmentation.csv").string()));
    ASSERT_NO_THROW(baysor::save_cell_stat_df(
        result, (output.path() / "segmentation_cell_stats.csv").string()));
    ASSERT_NO_THROW(baysor::save_matrix_to_tsv(
        result, (output.path() / "segmentation_counts.tsv").string()));
    ASSERT_NO_THROW(baysor::save_polygons_geojson(
        result,
        (output.path() / "segmentation_polygons_2d.json").string(),
        "FeatureCollection"));

    try {
        const auto labels = compare_molecule_tables(
            output.path() / "segmentation.csv",
            reference_dir / "segmentation.csv");
        compare_cell_statistics(
            output.path() / "segmentation_cell_stats.csv",
            reference_dir / "segmentation_cell_stats.csv",
            labels);
        compare_count_matrices(
            output.path() / "segmentation_counts.tsv",
            reference_dir / "segmentation_counts.tsv",
            labels);
        compare_polygons(
            output.path() / "segmentation_polygons_2d.json",
            reference_dir / "segmentation_polygons_2d.json",
            labels);
        compare_resolved_options(
            result.resolved_options,
            reference_dir / "resolved_params.toml");
    } catch (const std::exception& error) {
        FAIL() << error.what();
    }
}

struct ScientificParityCase {
    std::string name;
    fs::path molecules;
    fs::path config;
    fs::path reference;
    bool is_3d = false;
    bool has_clusters = false;
};

ScientificParityCase make_scientific_parity_case(
    std::string name,
    const fs::path& molecules,
    const fs::path& config,
    const fs::path& reference,
    bool is_3d,
    bool has_clusters
) {
    const fs::path fixture_root(BAYSOR_TEST_SCIENTIFIC_PARITY_DIR);
    return {
        std::move(name),
        fixture_root / molecules,
        fixture_root / config,
        fixture_root / reference,
        is_3d,
        has_clusters
    };
}

baysor::SegmentationRequest make_scientific_parity_request(
    const ScientificParityCase& test_case
) {
    const baysor::RunOptions config = baysor::load_config(test_case.config.string());
    baysor::SegmentationRequest request;
    request.molecules.path = test_case.molecules.string();
    request.molecules.options = config.molecules;
    request.prior = config.prior;
    request.segmentation = config.segmentation;
    request.neighborhood_composition.neighborhood_size =
        config.plotting.gene_composition_neighborhood;
    request.neighborhood_composition.method = config.plotting.ncv_method;
    request.requested_products = baysor::SegmentationProducts::none();
    request.requested_products.molecule_assignments = true;
    request.requested_products.molecule_confidence = true;
    request.requested_products.assignment_confidence = true;
    request.requested_products.molecule_clusters = test_case.has_clusters;
    request.requested_products.cell_statistics = true;
    request.requested_products.boundaries = true;
    request.requested_products.count_matrix = true;
    request.execution.native_threads = 1;
    return request;
}

void serialize_scientific_products(
    const baysor::SegmentationResult& result,
    const fs::path& output_dir
) {
    const auto paths = baysor::get_output_paths(
        output_dir.string(), baysor::OutputStyle::Legacy, "tsv");
    baysor::save_segmented_df(result, paths.segmented_df);
    baysor::save_cell_stat_df(result, paths.cell_stats);
    baysor::save_matrix_to_tsv(result, paths.counts);
    if (result.boundaries_3d) {
        baysor::save_polygon_stack_geojson(*result.boundaries_3d, paths, "FeatureCollection");
    } else {
        baysor::save_polygons_geojson(result, paths.polygons_2d, "FeatureCollection");
    }
}

void compare_scientific_artifacts(
    const fs::path& actual_dir,
    const fs::path& reference_dir,
    bool is_3d,
    bool has_clusters
) {
    CellLabelMap cluster_labels;
    const auto labels = compare_molecule_tables(
        actual_dir / "segmentation.csv",
        reference_dir / "segmentation.csv",
        has_clusters ? &cluster_labels : nullptr);
    compare_cell_statistics(
        actual_dir / "segmentation_cell_stats.csv",
        reference_dir / "segmentation_cell_stats.csv",
        labels,
        has_clusters ? &cluster_labels : nullptr);
    compare_count_matrices(
        actual_dir / "segmentation_counts.tsv",
        reference_dir / "segmentation_counts.tsv",
        labels);
    compare_polygons(
        actual_dir / "segmentation_polygons_2d.json",
        reference_dir / "segmentation_polygons_2d.json",
        labels);
    if (is_3d) {
        compare_polygon_stacks(
            actual_dir / "segmentation_polygons_3d.json",
            reference_dir / "segmentation_polygons_3d.json",
            labels);
    }
}

const baysor::SegmentationResult& require_completed(
    const baysor::SegmentationOutcome& outcome
) {
    require(std::holds_alternative<baysor::SegmentationResult>(outcome),
            "Segmentation was unexpectedly cancelled");
    return std::get<baysor::SegmentationResult>(outcome);
}

void compare_owned_polygons(
    const baysor::PolygonCollection& actual,
    const baysor::PolygonCollection& expected,
    const std::string& context
) {
    require(actual.size() == expected.size(), context + " polygon count changed");
    for (const auto& [cell, expected_vertices] : expected) {
        auto actual_it = actual.find(cell);
        require(actual_it != actual.end(), context + " is missing cell " + cell);
        require(actual_it->second.rows() == expected_vertices.rows() &&
                    actual_it->second.cols() == expected_vertices.cols(),
                context + " shape changed for cell " + cell);
        require(actual_it->second.isApprox(expected_vertices, 0.0),
                context + " coordinates changed for cell " + cell);
    }
}

void compare_owned_scientific_results(
    const baysor::SegmentationResult& actual,
    const baysor::SegmentationResult& expected
) {
    require(actual.molecules.source_transcript_id == expected.molecules.source_transcript_id,
            "Retained transcript identities changed");
    require(actual.cell_assignments == expected.cell_assignments, "Cell assignments changed");
    require(actual.molecule_clusters == expected.molecule_clusters, "Molecule clusters changed");
    require(actual.molecules.confidence == expected.molecules.confidence,
            "Molecule confidence changed");
    require(actual.assignment_confidence == expected.assignment_confidence,
            "Assignment confidence changed");
    require(actual.cell_ids == expected.cell_ids, "Cell identifiers changed");

    require(actual.cell_statistics.has_value() == expected.cell_statistics.has_value(),
            "Cell-statistics availability changed");
    if (actual.cell_statistics) {
        require(actual.cell_statistics->cell_ids == expected.cell_statistics->cell_ids,
                "Cell-statistics identifiers changed");
        require(actual.cell_statistics->columns == expected.cell_statistics->columns,
                "Cell-statistics columns changed");
        require(actual.cell_statistics->values.isApprox(expected.cell_statistics->values, 0.0),
                "Cell-statistics values changed");
    }

    require(actual.count_matrix.has_value() == expected.count_matrix.has_value(),
            "Count-matrix availability changed");
    if (actual.count_matrix) {
        require(actual.count_matrix->cell_ids == expected.count_matrix->cell_ids,
                "Count-matrix cell identifiers changed");
        require(actual.count_matrix->gene_names == expected.count_matrix->gene_names,
                "Count-matrix gene names changed");
        const Eigen::MatrixXf actual_counts(actual.count_matrix->values);
        const Eigen::MatrixXf expected_counts(expected.count_matrix->values);
        require(actual_counts.isApprox(expected_counts, 0.0), "Count-matrix values changed");
    }

    require(actual.boundaries_2d.has_value() == expected.boundaries_2d.has_value(),
            "2D-boundary availability changed");
    if (actual.boundaries_2d) {
        compare_owned_polygons(*actual.boundaries_2d, *expected.boundaries_2d, "2D boundaries");
    }
    require(actual.boundaries_3d.has_value() == expected.boundaries_3d.has_value(),
            "3D-boundary availability changed");
    if (actual.boundaries_3d) {
        require(actual.boundaries_3d->size() == expected.boundaries_3d->size(),
                "3D-boundary layer count changed");
        for (std::size_t layer = 0; layer < actual.boundaries_3d->size(); ++layer) {
            require(actual.boundaries_3d->at(layer).first == expected.boundaries_3d->at(layer).first,
                    "3D-boundary layer name changed");
            compare_owned_polygons(
                actual.boundaries_3d->at(layer).second,
                expected.boundaries_3d->at(layer).second,
                "3D boundaries in layer " + actual.boundaries_3d->at(layer).first);
        }
    }
}

// Fixture intent, semantic comparison rules, and historical-reference policy:
// tests/fixtures/scientific_parity/README.md
class ScientificParityTest : public ::testing::TestWithParam<ScientificParityCase> {};

TEST_P(ScientificParityTest, CurrentCliMatchesPreExtractionReference) {
    const auto& test_case = GetParam();
    TemporaryDirectory output;
    const std::string command =
        quote_argument(fs::path(BAYSOR_TEST_CMAKE_COMMAND)) +
        " -E env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE " +
        quote_argument(fs::path(BAYSOR_TEST_CLI_PATH)) +
        " run " + quote_argument(test_case.molecules) +
        " --config " + quote_argument(test_case.config) +
        " --output " + quote_argument(output.path()) +
        " --output-style legacy --polygon-format FeatureCollection" +
        " --count-matrix-format tsv --skip-ncv-color";

    ASSERT_EQ(std::system(command.c_str()), 0) << "Scientific-parity CLI command failed: " << command;
    try {
        compare_scientific_artifacts(
            output.path(), test_case.reference, test_case.is_3d, test_case.has_clusters);
        compare_resolved_options(
            output.path() / "segmentation_params.dump.toml",
            test_case.reference / "resolved_params.toml");
    } catch (const std::exception& error) {
        FAIL() << error.what();
    }
}

TEST_P(ScientificParityTest, DirectOperationMatchesPreExtractionReference) {
    const auto& test_case = GetParam();
    const auto outcome = baysor::run_segmentation(
        make_scientific_parity_request(test_case), baysor::CancellationToken{});
    ASSERT_TRUE(std::holds_alternative<baysor::SegmentationResult>(outcome));
    const auto& result = std::get<baysor::SegmentationResult>(outcome);
    EXPECT_EQ(result.molecules.is_3d(), test_case.is_3d);
    EXPECT_EQ(result.produced_products.molecule_clusters, test_case.has_clusters);

    TemporaryDirectory output;
    ASSERT_NO_THROW(serialize_scientific_products(result, output.path()));
    try {
        compare_scientific_artifacts(
            output.path(), test_case.reference, test_case.is_3d, test_case.has_clusters);
        compare_resolved_options(result.resolved_options, test_case.reference / "resolved_params.toml");
    } catch (const std::exception& error) {
        FAIL() << error.what();
    }
}

INSTANTIATE_TEST_SUITE_P(
    N2b,
    ScientificParityTest,
    ::testing::Values(
        make_scientific_parity_case(
            "DuplicateCoordinates2D",
            "duplicate_2d/molecules.csv",
            "duplicate_2d/config.toml",
            "duplicate_2d/reference",
            false,
            false),
        make_scientific_parity_case(
            "MrfClustering2D",
            "clustering_2d/molecules.csv",
            "clustering_2d/mrf.toml",
            "clustering_2d/reference/mrf",
            false,
            true),
        make_scientific_parity_case(
            "LouvainClustering2D",
            "clustering_2d/molecules.csv",
            "clustering_2d/louvain.toml",
            "clustering_2d/reference/louvain",
            false,
            true),
        make_scientific_parity_case(
            "LeidenClustering2D",
            "clustering_2d/molecules.csv",
            "clustering_2d/leiden.toml",
            "clustering_2d/reference/leiden",
            false,
            true),
        make_scientific_parity_case(
            "ThreeDimensional",
            "three_d/molecules.csv",
            "three_d/config.toml",
            "three_d/reference",
            true,
            false)),
    [](const ::testing::TestParamInfo<ScientificParityCase>& info) {
        return info.param.name;
    });

TEST(LegacySeedCompatibility, DefaultAndExplicitSeedUseTheSameScientificStream) {
    Eigen::MatrixXd duplicate_points(2, 4);
    duplicate_points << 0.0, 0.0, 1.0, 1.0,
                        0.0, 0.0, 1.0, 1.0;
    baysor::reset_global_xoshiro_rng(baysor::kDefaultSegmentationSeed);
    const Eigen::MatrixXd legacy_normalized = baysor::normalize_points(duplicate_points);
    baysor::Xoshiro256pp explicit_random_state(baysor::kDefaultSegmentationSeed);
    const Eigen::MatrixXd explicit_normalized =
        baysor::normalize_points(duplicate_points, &explicit_random_state);
    EXPECT_TRUE(legacy_normalized.isApprox(explicit_normalized, 0.0));
    EXPECT_EQ(baysor::global_xoshiro_rng()(), explicit_random_state());
    baysor::reset_global_xoshiro_rng(baysor::kDefaultSegmentationSeed);

    const auto test_case = make_scientific_parity_case(
        "DuplicateCoordinates2D",
        "duplicate_2d/molecules.csv",
        "duplicate_2d/config.toml",
        "duplicate_2d/reference",
        false,
        false);
    auto default_request = make_scientific_parity_request(test_case);
    auto explicit_request = default_request;
    explicit_request.random_seed = baysor::kDefaultSegmentationSeed;
    const auto default_outcome = baysor::run_segmentation(
        default_request, baysor::CancellationToken{});
    const auto explicit_outcome = baysor::run_segmentation(
        explicit_request, baysor::CancellationToken{});
    try {
        compare_owned_scientific_results(
            require_completed(default_outcome), require_completed(explicit_outcome));
    } catch (const std::exception& error) {
        FAIL() << error.what();
    }
}

bool is_well_formed_hex_color(const std::string& color) {
    return color.size() == 7 && color.front() == '#' &&
        std::all_of(color.begin() + 1, color.end(), [](unsigned char character) {
            return std::isxdigit(character) != 0;
        });
}

TEST(NeighborhoodCompositionColors, AreRepeatableAndScientificallyNonInterfering) {
    const auto test_case = make_scientific_parity_case(
        "LouvainClustering2D",
        "clustering_2d/molecules.csv",
        "clustering_2d/louvain.toml",
        "clustering_2d/reference/louvain",
        false,
        true);
    auto without_colors_request = make_scientific_parity_request(test_case);
    auto with_colors_request = without_colors_request;
    with_colors_request.requested_products.neighborhood_composition_colors = true;

    const auto without_colors_outcome = baysor::run_segmentation(
        without_colors_request, baysor::CancellationToken{});
    const auto first_colors_outcome = baysor::run_segmentation(
        with_colors_request, baysor::CancellationToken{});
    const auto second_colors_outcome = baysor::run_segmentation(
        with_colors_request, baysor::CancellationToken{});

    try {
        const auto& without_colors = require_completed(without_colors_outcome);
        const auto& first_colors = require_completed(first_colors_outcome);
        const auto& second_colors = require_completed(second_colors_outcome);
        require(without_colors.neighborhood_composition_colors.empty(),
                "An omitted NCV-colour product was unexpectedly materialized");
        require(first_colors.neighborhood_composition_colors.size() ==
                    static_cast<std::size_t>(first_colors.molecules.n_molecules()),
                "NCV colour count does not match the retained molecule count");
        require(std::all_of(
                    first_colors.neighborhood_composition_colors.begin(),
                    first_colors.neighborhood_composition_colors.end(),
                    is_well_formed_hex_color),
                "NCV colours are not well-formed #RRGGBB values");
        require(first_colors.neighborhood_composition_colors ==
                    second_colors.neighborhood_composition_colors,
                "NCV colours are not repeatable for the locked one-thread run");
        compare_owned_scientific_results(first_colors, without_colors);
        compare_owned_scientific_results(second_colors, without_colors);
    } catch (const std::exception& error) {
        FAIL() << error.what();
    }
}

}  // namespace
