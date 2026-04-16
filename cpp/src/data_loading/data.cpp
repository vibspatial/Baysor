#include "baysor/data_loading/data.h"
#include "baysor/utils/general.h"

#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <arrow/csv/api.h>
#include <arrow/io/api.h>
#include <arrow/result.h>
#include <arrow/status.h>
#include <parquet/arrow/reader.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace baysor {

// ============================================================================
// Helpers: extract columns from Arrow Table
// ============================================================================

namespace {

/// Get file extension, lowercased
std::string file_extension(const std::string& path) {
    auto dot = path.rfind('.');
    if (dot == std::string::npos) return "";
    std::string ext = path.substr(dot + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

/// Macro-like helper: check Arrow status/result
#define ARROW_CHECK_OK(expr)                                                   \
    do {                                                                        \
        auto _s = (expr);                                                       \
        if (!_s.ok())                                                           \
            throw std::runtime_error(std::string("Arrow error: ") +            \
                                     _s.ToString());                            \
    } while (0)

template <typename T>
T arrow_unwrap(arrow::Result<T>&& result) {
    if (!result.ok()) {
        throw std::runtime_error(std::string("Arrow error: ") + result.status().ToString());
    }
    return std::move(*result);
}

/// Read a CSV file into an Arrow Table
std::shared_ptr<arrow::Table> read_csv_arrow(const std::string& path) {
    auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path));
    auto read_options = arrow::csv::ReadOptions::Defaults();
    auto parse_options = arrow::csv::ParseOptions::Defaults();
    auto convert_options = arrow::csv::ConvertOptions::Defaults();
    // Let Arrow auto-detect types
    auto reader = arrow_unwrap(
        arrow::csv::TableReader::Make(arrow::io::default_io_context(),
                                       input, read_options, parse_options, convert_options));
    return arrow_unwrap(reader->Read());
}

/// Read a Parquet file into an Arrow Table
std::shared_ptr<arrow::Table> read_parquet_arrow(const std::string& path) {
    auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path));
    parquet::arrow::FileReaderBuilder builder;
    ARROW_CHECK_OK(builder.Open(input));
    auto reader = arrow_unwrap(builder.Build());
    std::shared_ptr<arrow::Table> table;
    ARROW_CHECK_OK(reader->ReadTable(&table));
    return table;
}

/// Read a tabular file (detect format by extension)
std::shared_ptr<arrow::Table> read_table(const std::string& path) {
    std::string ext = file_extension(path);
    if (ext == "csv" || ext == "tsv") {
        return read_csv_arrow(path);
    } else if (ext == "parquet" || ext == "pq") {
        return read_parquet_arrow(path);
    } else {
        throw std::runtime_error("Unsupported file format: ." + ext +
                                 ". Provide a CSV or Parquet file.");
    }
}

/// Check that a column exists in the table
int find_column_index(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
    auto schema = table->schema();
    int idx = schema->GetFieldIndex(name);
    if (idx < 0) {
        throw std::runtime_error("Column '" + name + "' not found in the data. "
                                 "Available columns: " +
                                 schema->ToString());
    }
    return idx;
}

/// Extract a double column from an Arrow table (handles chunked arrays, int/float types)
std::vector<double> extract_double_column(const std::shared_ptr<arrow::Table>& table,
                                           const std::string& name) {
    int idx = find_column_index(table, name);
    auto chunked = table->column(idx);
    int64_t n = chunked->length();
    std::vector<double> result(n);

    int64_t offset = 0;
    for (int c = 0; c < chunked->num_chunks(); ++c) {
        auto chunk = chunked->chunk(c);
        int64_t len = chunk->length();

        // Handle different numeric types
        if (chunk->type_id() == arrow::Type::DOUBLE) {
            auto arr = std::static_pointer_cast<arrow::DoubleArray>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = arr->Value(i);
            }
        } else if (chunk->type_id() == arrow::Type::FLOAT) {
            auto arr = std::static_pointer_cast<arrow::FloatArray>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = static_cast<double>(arr->Value(i));
            }
        } else if (chunk->type_id() == arrow::Type::INT64) {
            auto arr = std::static_pointer_cast<arrow::Int64Array>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = static_cast<double>(arr->Value(i));
            }
        } else if (chunk->type_id() == arrow::Type::INT32) {
            auto arr = std::static_pointer_cast<arrow::Int32Array>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = static_cast<double>(arr->Value(i));
            }
        } else if (chunk->type_id() == arrow::Type::INT16) {
            auto arr = std::static_pointer_cast<arrow::Int16Array>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = static_cast<double>(arr->Value(i));
            }
        } else {
            // Try casting to double
            auto maybe_casted = arrow::compute::Cast(*chunk, arrow::TypeHolder(arrow::float64()));
            if (!maybe_casted.ok()) {
                throw std::runtime_error("Cannot convert column '" + name +
                                         "' to double: " + maybe_casted.status().ToString());
            }
            auto casted = std::static_pointer_cast<arrow::DoubleArray>(*maybe_casted);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = casted->Value(i);
            }
        }
        offset += len;
    }
    return result;
}

/// Extract a string column from an Arrow table (handles string, large_string, dict types)
std::vector<std::string> extract_string_column(const std::shared_ptr<arrow::Table>& table,
                                                const std::string& name) {
    int idx = find_column_index(table, name);
    auto chunked = table->column(idx);
    int64_t n = chunked->length();
    std::vector<std::string> result(n);

    int64_t offset = 0;
    for (int c = 0; c < chunked->num_chunks(); ++c) {
        auto chunk = chunked->chunk(c);
        int64_t len = chunk->length();

        if (chunk->type_id() == arrow::Type::STRING) {
            auto arr = std::static_pointer_cast<arrow::StringArray>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = arr->GetString(i);
            }
        } else if (chunk->type_id() == arrow::Type::LARGE_STRING) {
            auto arr = std::static_pointer_cast<arrow::LargeStringArray>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = arr->GetString(i);
            }
        } else if (chunk->type_id() == arrow::Type::DICTIONARY) {
            // Dictionary-encoded strings (common in Parquet)
            auto dict_arr = std::static_pointer_cast<arrow::DictionaryArray>(chunk);
            auto dict_values = std::static_pointer_cast<arrow::StringArray>(dict_arr->dictionary());
            auto indices = dict_arr->indices();

            // Handle different index types
            if (indices->type_id() == arrow::Type::INT32) {
                auto idx_arr = std::static_pointer_cast<arrow::Int32Array>(indices);
                for (int64_t i = 0; i < len; ++i) {
                    if (!dict_arr->IsNull(i)) {
                        result[offset + i] = dict_values->GetString(idx_arr->Value(i));
                    }
                }
            } else if (indices->type_id() == arrow::Type::INT16) {
                auto idx_arr = std::static_pointer_cast<arrow::Int16Array>(indices);
                for (int64_t i = 0; i < len; ++i) {
                    if (!dict_arr->IsNull(i)) {
                        result[offset + i] = dict_values->GetString(idx_arr->Value(i));
                    }
                }
            } else if (indices->type_id() == arrow::Type::INT8) {
                auto idx_arr = std::static_pointer_cast<arrow::Int8Array>(indices);
                for (int64_t i = 0; i < len; ++i) {
                    if (!dict_arr->IsNull(i)) {
                        result[offset + i] = dict_values->GetString(idx_arr->Value(i));
                    }
                }
            } else {
                // Fallback: cast indices
                auto idx_arr = std::static_pointer_cast<arrow::Int64Array>(
                    arrow_unwrap(arrow::compute::Cast(*indices, arrow::TypeHolder(arrow::int64()))));
                for (int64_t i = 0; i < len; ++i) {
                    if (!dict_arr->IsNull(i)) {
                        result[offset + i] = dict_values->GetString(
                            static_cast<int>(idx_arr->Value(i)));
                    }
                }
            }
        } else {
            // Numeric column used as gene — convert to string
            for (int64_t i = 0; i < len; ++i) {
                auto scalar = arrow_unwrap(chunk->GetScalar(i));
                result[offset + i] = scalar->ToString();
            }
        }
        offset += len;
    }
    return result;
}

/// Check if a column exists
bool has_column(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
    return table->schema()->GetFieldIndex(name) >= 0;
}

} // anonymous namespace

// ============================================================================
// MoleculeData methods
// ============================================================================

Eigen::MatrixXd MoleculeData::position_matrix() const {
    int n = n_molecules();
    int d = n_dims();
    Eigen::MatrixXd mat(d, n);
    for (int i = 0; i < n; ++i) {
        mat(0, i) = x[i];
        mat(1, i) = y[i];
        if (d == 3) mat(2, i) = z[i];
    }
    return mat;
}

// ============================================================================
// Gene encoding
// ============================================================================

void encode_genes(MoleculeData& data, const std::vector<std::string>& gene_strings) {
    // Collect unique names and sort alphabetically
    std::set<std::string> unique_set;
    for (const auto& g : gene_strings) {
        if (!g.empty()) unique_set.insert(g);
    }

    data.gene_names.assign(unique_set.begin(), unique_set.end());

    // Build name -> 1-based ID map
    std::unordered_map<std::string, int> gene_id_map;
    gene_id_map.reserve(data.gene_names.size());
    for (int i = 0; i < static_cast<int>(data.gene_names.size()); ++i) {
        gene_id_map[data.gene_names[i]] = i + 1; // 1-based
    }

    // Encode
    data.gene.resize(gene_strings.size());
    for (size_t i = 0; i < gene_strings.size(); ++i) {
        auto it = gene_id_map.find(gene_strings[i]);
        data.gene[i] = (it != gene_id_map.end()) ? it->second : 0;
    }
}

// ============================================================================
// Gene filtering
// ============================================================================

void filter_genes_by_count(MoleculeData& data, int min_molecules_per_gene) {
    if (min_molecules_per_gene <= 1) return; // nothing to filter

    int n_genes = data.n_genes();
    int n = data.n_molecules();

    // Count molecules per gene
    std::vector<int> gene_counts(n_genes, 0);
    for (int i = 0; i < n; ++i) {
        int g = data.gene[i];
        if (g > 0 && g <= n_genes) gene_counts[g - 1]++;
    }

    // Find genes passing the threshold
    std::unordered_set<int> passing_genes;
    for (int g = 0; g < n_genes; ++g) {
        if (gene_counts[g] >= min_molecules_per_gene) {
            passing_genes.insert(g + 1); // 1-based
        }
    }

    if (static_cast<int>(passing_genes.size()) == n_genes) return; // all pass

    // Filter molecules: keep only those with passing genes
    std::vector<int> keep_indices;
    keep_indices.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (passing_genes.count(data.gene[i])) {
            keep_indices.push_back(i);
        }
    }

    // Compact data
    auto compact = [&](auto& vec) {
        if (vec.empty()) return;
        using T = typename std::decay<decltype(vec)>::type;
        T new_vec;
        new_vec.reserve(keep_indices.size());
        for (int idx : keep_indices) {
            new_vec.push_back(std::move(vec[idx]));
        }
        vec = std::move(new_vec);
    };

    compact(data.x);
    compact(data.y);
    compact(data.z);
    compact(data.confidence);
    compact(data.cluster);
    compact(data.prior_segmentation);
    compact(data.nuclei_probs);

    // Re-encode genes: collect raw strings for kept molecules, then re-encode
    std::vector<std::string> kept_gene_strings;
    kept_gene_strings.reserve(keep_indices.size());
    for (int idx : keep_indices) {
        int g = data.gene[idx];
        kept_gene_strings.push_back(data.gene_names[g - 1]);
    }

    // Clear old gene data and re-encode
    data.gene.clear();
    data.gene_names.clear();
    encode_genes(data, kept_gene_strings);
}

void filter_genes_by_pattern(MoleculeData& data, const std::vector<std::string>& patterns) {
    if (patterns.empty()) return;

    int n_genes = data.n_genes();
    int n = data.n_molecules();

    // Convert wildcard patterns to regexes
    // Julia uses Regex() directly, where users write "Blank*" meaning "Blank.*"
    // We support shell-glob style: * -> .*, ? -> .
    std::vector<std::regex> regexes;
    regexes.reserve(patterns.size());
    for (const auto& pat : patterns) {
        // Escape regex special chars except * and ?
        std::string re_str;
        for (char c : pat) {
            switch (c) {
                case '*': re_str += ".*"; break;
                case '?': re_str += "."; break;
                case '.': re_str += "\\."; break;
                case '+': re_str += "\\+"; break;
                case '(': re_str += "\\("; break;
                case ')': re_str += "\\)"; break;
                case '[': re_str += "\\["; break;
                case ']': re_str += "\\]"; break;
                case '{': re_str += "\\{"; break;
                case '}': re_str += "\\}"; break;
                case '^': re_str += "\\^"; break;
                case '$': re_str += "\\$"; break;
                case '|': re_str += "\\|"; break;
                case '\\': re_str += "\\\\"; break;
                default: re_str += c; break;
            }
        }
        regexes.emplace_back(re_str, std::regex::ECMAScript);
    }

    // Find gene names matching any pattern
    std::unordered_set<int> excluded_genes; // 1-based
    std::vector<std::string> excluded_names;
    for (int g = 0; g < n_genes; ++g) {
        const auto& name = data.gene_names[g];
        for (const auto& re : regexes) {
            if (std::regex_match(name, re)) {
                excluded_genes.insert(g + 1);
                excluded_names.push_back(name);
                break;
            }
        }
    }

    if (excluded_genes.empty()) return;

    // Filter molecules: keep those NOT matching excluded genes
    std::vector<int> keep_indices;
    keep_indices.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (!excluded_genes.count(data.gene[i])) {
            keep_indices.push_back(i);
        }
    }

    // Compact (same pattern as filter_genes_by_count)
    auto compact = [&](auto& vec) {
        if (vec.empty()) return;
        using T = typename std::decay<decltype(vec)>::type;
        T new_vec;
        new_vec.reserve(keep_indices.size());
        for (int idx : keep_indices) {
            new_vec.push_back(std::move(vec[idx]));
        }
        vec = std::move(new_vec);
    };

    compact(data.x);
    compact(data.y);
    compact(data.z);
    compact(data.confidence);
    compact(data.cluster);
    compact(data.prior_segmentation);
    compact(data.nuclei_probs);

    // Re-encode
    std::vector<std::string> kept_gene_strings;
    kept_gene_strings.reserve(keep_indices.size());
    for (int idx : keep_indices) {
        int g = data.gene[idx];
        kept_gene_strings.push_back(data.gene_names[g - 1]);
    }

    data.gene.clear();
    data.gene_names.clear();
    encode_genes(data, kept_gene_strings);
}

// ============================================================================
// Read a single string column (for prior segmentation column references)
// ============================================================================

std::vector<std::string> read_string_column(const std::string& path, const std::string& col_name) {
    auto table = read_table(path);
    return extract_string_column(table, col_name);
}

// ============================================================================
// Main loading function
// ============================================================================

RawTableData read_tabular_file(const std::string& path, const DataOptions& opts) {
    auto table = read_table(path);
    RawTableData raw;

    // Extract x, y (required)
    raw.x = extract_double_column(table, opts.x_col);
    raw.y = extract_double_column(table, opts.y_col);

    // Extract z (optional)
    if (!opts.force_2d && has_column(table, opts.z_col)) {
        raw.z = extract_double_column(table, opts.z_col);
        raw.has_z = true;

        // Drop z if only 1 unique value
        if (!raw.z.empty()) {
            double z0 = raw.z[0];
            bool all_same = true;
            for (size_t i = 1; i < raw.z.size(); ++i) {
                if (raw.z[i] != z0) { all_same = false; break; }
            }
            if (all_same) {
                raw.z.clear();
                raw.has_z = false;
            }
        }
    }

    // Extract gene column as strings
    raw.gene_str = extract_string_column(table, opts.gene_col);

    // Stringify any non-string gene values (handles numeric gene columns)
    // Already handled in extract_string_column via scalar->ToString()

    return raw;
}

MoleculeData load_molecules(const std::string& path, const DataOptions& opts) {
    MoleculeData data;

    // Read the full table once (shared with optional column extraction below)
    auto table = read_table(path);

    // Extract required spatial columns
    data.x = extract_double_column(table, opts.x_col);
    data.y = extract_double_column(table, opts.y_col);

    // Extract z (optional)
    if (!opts.force_2d && has_column(table, opts.z_col)) {
        data.z = extract_double_column(table, opts.z_col);
        // Drop z if only 1 unique value
        if (!data.z.empty()) {
            double z0 = data.z[0];
            bool all_same = true;
            for (size_t i = 1; i < data.z.size(); ++i) {
                if (data.z[i] != z0) { all_same = false; break; }
            }
            if (all_same) {
                data.z.clear();
            }
        }
    }

    // Extract gene column and encode
    auto gene_strings = extract_string_column(table, opts.gene_col);
    encode_genes(data, gene_strings);

    // Load optional metadata columns (if present in the input file)
    if (has_column(table, "confidence")) {
        data.confidence = extract_double_column(table, "confidence");
    }
    if (has_column(table, "cluster")) {
        // cluster column may be int or string; read as double and round
        auto cluster_dbl = extract_double_column(table, "cluster");
        data.cluster.resize(cluster_dbl.size());
        for (size_t i = 0; i < cluster_dbl.size(); ++i) {
            data.cluster[i] = static_cast<int>(std::round(cluster_dbl[i]));
        }
    }
    if (has_column(table, "nuclei_probs")) {
        data.nuclei_probs = extract_double_column(table, "nuclei_probs");
    }

    // Filter genes with too few molecules
    if (opts.min_molecules_per_gene > 1) {
        filter_genes_by_count(data, opts.min_molecules_per_gene);
    }

    // Filter genes matching exclude patterns
    if (!opts.exclude_genes.empty()) {
        auto patterns = split_string_list(opts.exclude_genes);
        if (!patterns.empty()) {
            filter_genes_by_pattern(data, patterns);
        }
    }

    return data;
}

} // namespace baysor
