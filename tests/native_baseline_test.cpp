#include <gtest/gtest.h>

#include "baysor/utils/options.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <chrono>
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
    std::size_t consumed = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid number '" + value + "' for " + context);
    }
    require(consumed == value.size(), "Invalid number '" + value + "' for " + context);
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

CellLabelMap compare_molecule_tables(const fs::path& actual_path, const fs::path& reference_path) {
    const auto actual = read_delimited_table(actual_path, ',');
    const auto reference = read_delimited_table(reference_path, ',');
    require(actual.columns == reference.columns, "Segmented-molecule columns changed");

    const auto actual_rows = rows_by_key(actual, "transcript_id");
    const auto reference_rows = rows_by_key(reference, "transcript_id");
    require(actual_rows.size() == reference_rows.size(), "Retained transcript count changed");

    CellLabelMap labels;
    std::map<std::string, std::string> reverse_labels;
    for (const auto& [transcript_id, reference_row] : reference_rows) {
        auto actual_it = actual_rows.find(transcript_id);
        require(actual_it != actual_rows.end(), "Missing transcript_id " + transcript_id);
        const auto& actual_row = *actual_it->second;

        require(actual.value(actual_row, "gene") == reference.value(*reference_row, "gene"),
                "Gene changed for transcript_id " + transcript_id);
        for (const auto& coordinate : {"x", "y"}) {
            require_near(
                parse_number(actual.value(actual_row, coordinate), coordinate),
                parse_number(reference.value(*reference_row, coordinate), coordinate),
                1e-9,
                1e-9,
                std::string(coordinate) + " for transcript_id " + transcript_id
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
            auto [mapping, inserted] = labels.emplace(actual_cell, reference_cell);
            require(inserted || mapping->second == reference_cell,
                    "Candidate cell " + actual_cell + " combines reference partitions");
            auto [reverse, reverse_inserted] = reverse_labels.emplace(reference_cell, actual_cell);
            require(reverse_inserted || reverse->second == actual_cell,
                    "Reference cell " + reference_cell + " was split across candidate cells");
        }
    }
    require(labels.size() == reverse_labels.size(), "Cell relabelling is not bijective");
    return labels;
}

void compare_cell_statistics(
    const fs::path& actual_path,
    const fs::path& reference_path,
    const CellLabelMap& labels
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

        for (const auto& column : {
                 "x", "y", "n_transcripts", "density", "elongation", "area",
                 "avg_confidence", "avg_assignment_confidence", "lifespan"
             }) {
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

PolygonMap read_polygons(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open polygons: " + path.string());
    nlohmann::json document;
    input >> document;
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

void compare_polygons(
    const fs::path& actual_path,
    const fs::path& reference_path,
    const CellLabelMap& labels
) {
    const auto actual = read_polygons(actual_path);
    const auto reference = read_polygons(reference_path);
    require(actual.size() == reference.size(), "Boundary polygon count changed");

    for (const auto& [actual_cell, actual_vertices] : actual) {
        auto label_it = labels.find(actual_cell);
        require(label_it != labels.end(), "Boundaries contain unknown cell " + actual_cell);
        auto reference_it = reference.find(label_it->second);
        require(reference_it != reference.end(), "Reference boundaries are missing a cell");
        const auto& reference_vertices = reference_it->second;
        require(actual_vertices.size() == reference_vertices.size(),
                "Boundary vertex count changed for cell " + actual_cell);
        for (std::size_t i = 0; i < actual_vertices.size(); ++i) {
            for (std::size_t dimension = 0; dimension < 2; ++dimension) {
                require_near(
                    actual_vertices[i][dimension],
                    reference_vertices[i][dimension],
                    1e-8,
                    1e-8,
                    "boundary coordinate for cell " + actual_cell
                );
            }
        }
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

}  // namespace
