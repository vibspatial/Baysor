#include "baysor/reporting/output.h"
#include "baysor/data_loading/data.h"

#include <hdf5.h>

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

namespace baysor {

// ============================================================================
// save_matrix_to_loom  —  C API implementation (avoids C++ ABI issues)
// ============================================================================
//
// Loom spec: https://linnarssonlab.org/loompy/format/index.html
//
// Layout produced (mirrors Julia's save_matrix_to_loom):
//   /matrix            float32  n_cells × n_genes  (chunked, shuffle+deflate)
//   /row_attrs/Name    variable-length UTF-8 strings  [n_genes]
//   /col_attrs/Name    variable-length UTF-8 strings  [n_cells]
//   /col_attrs/CellID  float64  [n_cells]   (1-based indices)
//   /col_attrs/<key>   float64[] or vlen-string[] per extra attribute
//   /attrs/LOOM_SPEC_VERSION  "3.0.0"

// RAII wrappers so we don't leak handles on exceptions.
struct HidGuard {
    hid_t id;
    explicit HidGuard(hid_t h) : id(h) {}
    ~HidGuard() { if (id >= 0) H5Idec_ref(id); }
    operator hid_t() const { return id; }
};

static void check(hid_t id, const char* ctx) {
    if (id < 0) throw std::runtime_error(std::string("HDF5 error: ") + ctx);
}
static void check(herr_t err, const char* ctx) {
    if (err < 0) throw std::runtime_error(std::string("HDF5 error: ") + ctx);
}

// Write a 1-D array of variable-length UTF-8 strings.
static void write_vlen_strings(hid_t grp, const char* name,
                                const std::vector<std::string>& strs) {
    hid_t stype = H5Tcopy(H5T_C_S1);
    H5Tset_size(stype, H5T_VARIABLE);
    H5Tset_cset(stype, H5T_CSET_UTF8);

    hsize_t dim = strs.size();
    hid_t space = H5Screate_simple(1, &dim, nullptr);
    hid_t ds = H5Dcreate2(grp, name, stype, space,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    check(ds, name);

    std::vector<const char*> ptrs(strs.size());
    for (size_t i = 0; i < strs.size(); ++i) ptrs[i] = strs[i].c_str();
    check(H5Dwrite(ds, stype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrs.data()), name);

    H5Dclose(ds);
    H5Sclose(space);
    H5Tclose(stype);
}

// Write a 1-D array of float64.
static void write_doubles(hid_t grp, const char* name,
                           const std::vector<double>& vals) {
    hsize_t dim = vals.size();
    hid_t space = H5Screate_simple(1, &dim, nullptr);
    hid_t ds = H5Dcreate2(grp, name, H5T_NATIVE_DOUBLE, space,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    check(ds, name);
    check(H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals.data()), name);
    H5Dclose(ds);
    H5Sclose(space);
}

void save_matrix_to_loom(
    const Eigen::SparseMatrix<float>& matrix,   // n_cells × n_genes
    const std::vector<std::string>& gene_names,
    const std::vector<std::string>& cell_names,
    const std::string& path,
    const LoomColAttrs& col_attrs
) {
    int n_cells = static_cast<int>(matrix.rows());
    int n_genes = static_cast<int>(matrix.cols());

    if (static_cast<int>(gene_names.size()) != n_genes)
        throw std::runtime_error("save_matrix_to_loom: gene_names length mismatch");
    if (static_cast<int>(cell_names.size()) != n_cells)
        throw std::runtime_error("save_matrix_to_loom: cell_names length mismatch");

    // Open file
    hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    check(fid, "H5Fcreate");

    // ---- /matrix (float32, chunked, shuffle + deflate-3) ----
    {
        hsize_t dims[2]  = { static_cast<hsize_t>(n_cells),
                             static_cast<hsize_t>(n_genes) };
        hsize_t chunk[2] = { static_cast<hsize_t>(std::min(n_cells, 64)),
                             static_cast<hsize_t>(std::min(n_genes, 64)) };

        hid_t space  = H5Screate_simple(2, dims, nullptr);
        hid_t plist  = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist, 2, chunk);
        H5Pset_shuffle(plist);
        H5Pset_deflate(plist, 3);

        hid_t ds = H5Dcreate2(fid, "matrix", H5T_NATIVE_FLOAT, space,
                               H5P_DEFAULT, plist, H5P_DEFAULT);
        check(ds, "create /matrix");

        // Write row-by-row (one dense cell vector at a time) to avoid
        // materialising the full dense matrix.
        Eigen::SparseMatrix<float, Eigen::RowMajor> rm(matrix);
        std::vector<float> row_buf(n_genes, 0.0f);

        hsize_t row_dims[2]   = {1, static_cast<hsize_t>(n_genes)};
        hsize_t row_offset[2] = {0, 0};
        hid_t mem_space = H5Screate_simple(2, row_dims, nullptr);

        for (int i = 0; i < n_cells; ++i) {
            std::fill(row_buf.begin(), row_buf.end(), 0.0f);
            for (Eigen::SparseMatrix<float, Eigen::RowMajor>::InnerIterator it(rm, i);
                 it; ++it) {
                row_buf[it.col()] = it.value();
            }
            row_offset[0] = static_cast<hsize_t>(i);
            hid_t file_space = H5Dget_space(ds);
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET,
                                row_offset, nullptr, row_dims, nullptr);
            H5Dwrite(ds, H5T_NATIVE_FLOAT, mem_space, file_space, H5P_DEFAULT,
                     row_buf.data());
            H5Sclose(file_space);
        }

        H5Sclose(mem_space);
        H5Dclose(ds);
        H5Pclose(plist);
        H5Sclose(space);
    }

    // ---- /attrs ----
    {
        hid_t grp = H5Gcreate2(fid, "attrs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(grp, "create /attrs");
        write_vlen_strings(grp, "LOOM_SPEC_VERSION", {"3.0.0"});
        H5Gclose(grp);
    }

    // ---- /row_attrs ----
    {
        hid_t grp = H5Gcreate2(fid, "row_attrs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(grp, "create /row_attrs");
        write_vlen_strings(grp, "Name", gene_names);
        H5Gclose(grp);
    }

    // ---- /col_attrs ----
    {
        hid_t grp = H5Gcreate2(fid, "col_attrs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(grp, "create /col_attrs");

        write_vlen_strings(grp, "Name", cell_names);

        std::vector<double> cell_ids(n_cells);
        for (int i = 0; i < n_cells; ++i) cell_ids[i] = i + 1.0;
        write_doubles(grp, "CellID", cell_ids);

        for (const auto& [key, val] : col_attrs) {
            if (std::holds_alternative<std::vector<std::string>>(val)) {
                write_vlen_strings(grp, key.c_str(),
                                   std::get<std::vector<std::string>>(val));
            } else {
                write_doubles(grp, key.c_str(),
                              std::get<std::vector<double>>(val));
            }
        }

        H5Gclose(grp);
    }

    H5Fclose(fid);
}

// ============================================================================
// get_output_paths
// ============================================================================

OutputPaths get_output_paths(const std::string& base_path,
                              const std::string& count_matrix_format) {
    // Strip trailing slash if present
    std::string base = base_path;
    while (!base.empty() && (base.back() == '/' || base.back() == '\\')) base.pop_back();

    OutputPaths p;
    p.segmented_df       = base + "/segmentation.csv";
    p.cell_stats         = base + "/segmentation_cell_stats.csv";
    p.diagnostic_report  = base + "/diagnostic_report.html";
    p.molecule_plot      = base + "/segmentation_plot.html";
    p.polygons_2d        = base + "/segmentation_polygons_2d.json";
    p.polygons_3d        = base + "/segmentation_polygons_3d.json";
    p.params_dump        = base + "/segmentation_params.dump.toml";
    p.log_file           = base + "/segmentation_log.log";

    if (count_matrix_format == "tsv") {
        p.counts = base + "/segmentation_counts.tsv";
    } else {
        p.counts = base + "/segmentation_counts.loom";
    }
    return p;
}

// ============================================================================
// save_segmented_df — per-molecule CSV
// ============================================================================

void save_segmented_df(const MoleculeData& data,
                       const std::vector<int>& assignment,
                       const std::vector<std::string>& gene_names,
                       const std::string& path,
                       const std::vector<std::string>* ncv_color,
                       const std::vector<double>* assignment_confidence,
                       const std::vector<int>* cluster) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_segmented_df: cannot open " + path);

    bool has_z    = data.is_3d();
    bool has_conf = !data.confidence.empty();
    bool has_col  = ncv_color && !ncv_color->empty();
    bool has_ac   = assignment_confidence && !assignment_confidence->empty();
    bool has_cl   = cluster && !cluster->empty();

    // Header
    f << "cell,gene,x,y";
    if (has_z)    f << ",z";
    if (has_conf) f << ",confidence";
    if (has_cl)   f << ",cluster";
    if (has_col)  f << ",ncv_color";
    if (has_ac)   f << ",assignment_confidence";
    f << ",is_noise\n";

    int n = data.n_molecules();
    for (int i = 0; i < n; ++i) {
        int cell = assignment[i];
        std::string cell_name = (cell > 0) ? ("cell_" + std::to_string(cell)) : "0";
        int g = data.gene[i];
        std::string gene_name = (g > 0 && g <= static_cast<int>(gene_names.size()))
                                ? gene_names[g - 1] : std::to_string(g);

        f << cell_name << ',' << gene_name << ',' << data.x[i] << ',' << data.y[i];
        if (has_z)    f << ',' << data.z[i];
        if (has_conf) f << ',' << data.confidence[i];
        if (has_cl)   f << ',' << (*cluster)[i];
        if (has_col)  f << ',' << (*ncv_color)[i];
        if (has_ac)   f << ',' << (*assignment_confidence)[i];
        f << ',' << (cell == 0 ? 1 : 0) << '\n';
    }
}

// ============================================================================
// save_cell_stat_df — per-cell CSV
// ============================================================================

void save_cell_stat_df(const Eigen::MatrixXd& stats,
                       const std::vector<std::string>& cell_names,
                       const std::vector<std::string>& col_names,
                       const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_cell_stat_df: cannot open " + path);

    f << "cell";
    for (const auto& c : col_names) f << ',' << c;
    f << '\n';

    int nr = static_cast<int>(stats.rows());
    int nc = static_cast<int>(stats.cols());
    for (int r = 0; r < nr; ++r) {
        f << cell_names[r];
        for (int c = 0; c < nc; ++c) f << ',' << stats(r, c);
        f << '\n';
    }
}

// ============================================================================
// save_matrix_to_tsv
// ============================================================================

void save_matrix_to_tsv(const Eigen::SparseMatrix<double>& matrix,
                         const std::vector<std::string>& gene_names,
                         const std::vector<std::string>& cell_names,
                         const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_matrix_to_tsv: cannot open " + path);

    // Header: cell names
    f << "gene";
    for (const auto& cn : cell_names) f << '\t' << cn;
    f << '\n';

    // Convert to column-major for efficient row iteration
    Eigen::SparseMatrix<double, Eigen::RowMajor> rm(matrix);
    int n_genes = static_cast<int>(rm.rows());
    for (int gi = 0; gi < n_genes; ++gi) {
        f << gene_names[gi];
        // Dense row output (sparse matrix — iterate nonzeros + fill zeros)
        std::vector<double> row(cell_names.size(), 0.0);
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(rm, gi); it; ++it) {
            row[it.col()] = it.value();
        }
        for (double v : row) f << '\t' << v;
        f << '\n';
    }
}

// ============================================================================
// save_polygons_geojson
// ============================================================================

void save_polygons_geojson(const PolygonCollection& polygons,
                            const std::string& path,
                            const std::string& format) {
    if (format == "none" || polygons.empty()) return;

    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_polygons_geojson: cannot open " + path);

    // Determine wrapper type (FeatureCollection or GeometryCollection)
    bool is_feature = (format != "GeometryCollection");

    if (is_feature) {
        f << "{\"type\":\"FeatureCollection\",\"features\":[\n";
    } else {
        f << "{\"type\":\"GeometryCollection\",\"geometries\":[\n";
    }

    bool first = true;
    for (const auto& [cell_name, hull] : polygons) {
        if (!first) f << ",\n";
        first = false;

        // Hull is 2×n (row 0 = x, row 1 = y); may be empty for tiny cells
        int nv = (hull.cols() > 0) ? static_cast<int>(hull.cols()) : 0;

        // Build coordinate ring (GeoJSON rings must be closed: first == last)
        std::string coords = "[[";
        if (nv > 0) {
            for (int j = 0; j < nv; ++j) {
                if (j > 0) coords += ",";
                coords += "[" + std::to_string(hull(0, j)) + ","
                              + std::to_string(hull(1, j)) + "]";
            }
            // Close the ring by repeating the first vertex
            coords += ",[" + std::to_string(hull(0, 0)) + ","
                           + std::to_string(hull(1, 0)) + "]";
        }
        coords += "]]";

        std::string geom = "{\"type\":\"Polygon\",\"coordinates\":" + coords + "}";

        if (is_feature) {
            f << "{\"type\":\"Feature\","
              << "\"geometry\":" << geom << ","
              << "\"properties\":{\"cell\":\"" << cell_name << "\"}}";
        } else {
            f << geom;
        }
    }

    if (is_feature) {
        f << "\n]}";
    } else {
        f << "\n]}";
    }
}

} // namespace baysor
