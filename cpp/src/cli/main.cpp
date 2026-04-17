#include <third_party/CLI11.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "baysor/utils/options.h"
#include "baysor/data_loading/data.h"
#include "baysor/data_loading/prior_segmentation.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/processing/data_processing/initialization.h"
#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/processing/utils/convex_hull.h"
#include "baysor/processing/bmm_algorithm/tracing.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/reporting/output.h"
#include "baysor/reporting/preview_report.h"

#include "baysor/utils/general.h"

#include <Eigen/Dense>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>

using namespace baysor;

// ============================================================================
// Subcommand: run
// ============================================================================

int cmd_run(
    const std::string& coordinates,
    const std::string& prior_segmentation,
    RunOptions& opts,
    const std::string& output,
    bool plot,
    bool skip_ncv_color,
    const std::string& polygon_format,
    const std::string& count_matrix_format,
    const std::string& cli_cmd
) {
    spdlog::info("Loading data from '{}'...", coordinates);

    fill_and_check_data_options(opts.data);

    auto data = load_molecules(coordinates, opts.data);
    spdlog::info("Loaded {} transcripts, {} genes.", data.n_molecules(), data.n_genes());

    fill_and_check_plotting_options(opts.plotting, opts.data.min_molecules_per_cell, data.n_genes());

    if (opts.segmentation.n_cells_init <= 0) {
        opts.segmentation.n_cells_init = default_param_value(
            "n_cells_init", opts.data.min_molecules_per_cell, data.n_molecules());
    }

    // Load prior segmentation if provided
    if (!prior_segmentation.empty()) {
        auto [scale, scale_std] = load_prior_segmentation(
            data, prior_segmentation, coordinates,
            opts.segmentation.unassigned_prior_label,
            opts.data.min_molecules_per_segment,
            opts.data.min_molecules_per_cell,
            opts.segmentation.estimate_scale_from_centers
        );

        if (opts.segmentation.estimate_scale_from_centers && scale > 0) {
            opts.segmentation.scale = scale;
            opts.segmentation.scale_std = std::to_string(scale_std);
        }
    }

    if (opts.segmentation.scale <= 0) {
        spdlog::error("Scale could not be determined. Either provide prior_segmentation or set --scale.");
        return 1;
    }

    spdlog::info("Using scale={:.2f}, scale_std={}",
                 opts.segmentation.scale, opts.segmentation.scale_std);

    double psc = opts.segmentation.prior_segmentation_confidence;

    spdlog::info("Estimating confidence...");
    append_confidence(data, opts.data.confidence_nn_id, psc);

    // Build molecule adjacency graph (MRF)
    spdlog::info("Building molecule graph...");
    auto adj_list = build_molecule_graph(data);

    // Create output directory
    {
        std::string mkdir_cmd = "mkdir -p \"" + output + "\"";
        if (std::system(mkdir_cmd.c_str()) != 0) {
            spdlog::warn("Could not create output directory '{}'", output);
        }
    }
    auto out_paths = get_output_paths(output, count_matrix_format);

    // Set up dual logger: console + log file (matches Julia's setup_logger)
    {
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        auto file_sink    = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
                                out_paths.log_file, /*truncate=*/true);
        auto logger = std::make_shared<spdlog::logger>(
                          "baysor", spdlog::sinks_init_list{console_sink, file_sink});
        logger->set_level(spdlog::level::info);
        spdlog::set_default_logger(logger);
    }

    // Guard unimplemented features
    if (!opts.segmentation.nuclei_genes.empty() || !opts.segmentation.cyto_genes.empty()) {
        spdlog::error("--nuclei-genes / --cyto-genes compartment segmentation is not yet implemented.");
        return 1;
    }

    int min_mols     = opts.data.min_molecules_per_cell;
    int n_iters   = opts.segmentation.iters;
    double scale  = opts.segmentation.scale;
    int n_cells   = opts.segmentation.n_cells_init;
    const std::string& scale_std = opts.segmentation.scale_std;

    // Optional molecule clustering (pre-segmentation cell type assignment)
    // Uses ICA initialization (matching Julia), falls back to hash init on failure.
    std::vector<int> mol_clusters;
    if (opts.segmentation.n_clusters > 1) {
        spdlog::info("Clustering molecules into {} types (ICA init)...",
                     opts.segmentation.n_clusters);
        auto cl = cluster_molecules_ica(
            data.gene, adj_list, data.confidence,
            opts.segmentation.n_clusters,
            /*tol=*/0.01, /*mrf_weight=*/1.0, /*max_iters=*/-1, /*verbose=*/true);
        mol_clusters = std::move(cl.assignment);  // 1-based cluster IDs
        spdlog::info("Molecule clustering complete.");
    }

    // Dispatch on dimensionality
    auto run_segmentation = [&](auto tag) {
        constexpr int N = decltype(tag)::value;

        spdlog::info("Initializing BmmData ({}D)...", N);
        auto bm_data = initialize_bmm_data<N>(
            data, adj_list, n_cells, scale, scale_std, psc, min_mols, /*verbose=*/true);

        // Wire molecule clusters into BmmData
        if (!mol_clusters.empty()) {
            bm_data.cluster_per_molecule = mol_clusters;
        }

        // History depth: match Julia's round(iters * 0.1)
        int history_depth = std::max(1, n_iters / 10);

        spdlog::info("Running segmentation ({} iters, history_depth={}, tol={})...",
                     n_iters, history_depth, opts.segmentation.tol);
        // Julia hardcodes min_n_samples=2 in drop_unused_components! — match that exactly.
        // min_mols = min_molecules_per_cell = display threshold only.
        bmm(bm_data, /*min_molecules_drop=*/2, n_iters,
            history_depth,
            /*verbose=*/true,
            /*component_split_step=*/3,
            /*refine=*/true,
            /*freeze_composition=*/false,
            /*freeze_position=*/false,
            /*freeze_components=*/false,
            opts.segmentation.tol,
            /*min_molecules_display=*/min_mols);

        int n_cells_final = bm_data.n_components();
        spdlog::info("Segmentation complete: {} cells.", n_cells_final);

        std::vector<std::string> ncv_color;
        if (!skip_ncv_color) {
            spdlog::info("Computing neighborhood composition colors...");
            int comp_k = opts.plotting.gene_composition_neighborhood;
            auto pos = data.position_matrix();
            auto neighb_cm = neighborhood_count_matrix(pos, data.gene, comp_k, data.n_genes());
            for (int k = 0; k < neighb_cm.outerSize(); ++k)
                for (Eigen::SparseMatrix<float>::InnerIterator it(neighb_cm, k); it; ++it)
                    it.valueRef() = static_cast<float>(std::log(it.value() * 10000.0f + 1e-5f));
            auto mol_vecs = estimate_gene_vectors(neighb_cm, data.gene, 20, "ri", true);
            ncv_color = gene_composition_color_embedding(mol_vecs, data.confidence);
        }

        // Save per-molecule CSV
        spdlog::info("Saving segmented molecule table...");
        const std::vector<double>* ac_ptr = bm_data.assignment_confidence.empty()
                                            ? nullptr : &bm_data.assignment_confidence;
        save_segmented_df(data, bm_data.assignment, data.gene_names, out_paths.segmented_df,
                          &ncv_color, ac_ptr,
                          mol_clusters.empty() ? nullptr : &mol_clusters);

        // Save per-cell stats CSV
        spdlog::info("Saving cell stats...");
        {
            std::vector<std::string> cell_names(n_cells_final);
            for (int i = 0; i < n_cells_final; ++i)
                cell_names[i] = "cell_" + std::to_string(i + 1);

            auto ids_by_cell = split_ids(bm_data.assignment, n_cells_final, true);

            // Precompute lifespan map (guid → lifespan)
            std::unordered_map<int,int> lifespan_map;
            if (!bm_data.assignment_history.empty()) {
                lifespan_map = estimate_component_lifespan(bm_data.assignment_history);
            }

            // Column names — order matches Julia:
            // cell, x, y, [z,] [cluster,] n_transcripts, density, elongation, area,
            // avg_confidence, [avg_assignment_confidence,] [max_cluster_frac,] [lifespan]
            bool has_cluster = !bm_data.cluster_per_cell.empty();
            bool has_ac      = !bm_data.assignment_confidence.empty();
            bool has_clmol   = !bm_data.cluster_per_molecule.empty();
            bool has_lspan   = !lifespan_map.empty();

            std::vector<std::string> col_names = {"x", "y"};
            if (N == 3) col_names.push_back("z");
            if (has_cluster) col_names.push_back("cluster");    // matches Julia position
            col_names.push_back("n_transcripts");
            col_names.push_back("density");
            col_names.push_back("elongation");
            col_names.push_back("area");
            col_names.push_back("avg_confidence");
            if (has_ac)    col_names.push_back("avg_assignment_confidence");
            if (has_clmol) col_names.push_back("max_cluster_frac");
            if (has_lspan) col_names.push_back("lifespan");

            // Build name→index map so writing order is independent of names order
            std::unordered_map<std::string, int> ci_map;
            for (int i = 0; i < static_cast<int>(col_names.size()); ++i)
                ci_map[col_names[i]] = i;

            int n_cols = static_cast<int>(col_names.size());
            Eigen::MatrixXd stats(n_cells_final, n_cols);
            stats.fill(std::numeric_limits<double>::quiet_NaN());

            for (int ci = 0; ci < n_cells_final; ++ci) {
                const auto& ids = ids_by_cell[ci];
                int np = static_cast<int>(ids.size());

                // --- Position means ---
                double sx = 0, sy = 0, sz = 0;
                double sum_conf = 0.0, sum_ac = 0.0;
                for (int mol : ids) {
                    sx += data.x[mol]; sy += data.y[mol];
                    if (N == 3 && !data.z.empty()) sz += data.z[mol];
                    sum_conf += data.confidence[mol];
                    if (has_ac) sum_ac += bm_data.assignment_confidence[mol];
                }
                double denom = np > 0 ? np : 1.0;

                stats(ci, ci_map["x"]) = sx / denom;
                stats(ci, ci_map["y"]) = sy / denom;
                if (N == 3) stats(ci, ci_map["z"]) = sz / denom;

                stats(ci, ci_map["n_transcripts"]) = np;
                stats(ci, ci_map["avg_confidence"]) = sum_conf / denom;
                if (has_ac)    stats(ci, ci_map["avg_assignment_confidence"]) = sum_ac / denom;
                if (has_cluster) stats(ci, ci_map["cluster"]) = bm_data.cluster_per_cell[ci];

                // --- Convex hull metrics (only for cells with > 2 molecules) ---
                if (np > 2) {
                    Eigen::MatrixXd pos2d(2, np);
                    for (int j = 0; j < np; ++j) {
                        pos2d(0, j) = data.x[ids[j]];
                        pos2d(1, j) = data.y[ids[j]];
                    }
                    auto hull = convex_hull(pos2d);
                    double area = polygon_area(hull);
                    stats(ci, ci_map["area"])    = area;
                    stats(ci, ci_map["density"]) = (area > 0) ? (np / area)
                                                               : std::numeric_limits<double>::quiet_NaN();

                    // Elongation: eigenvalue ratio of 2D sample covariance
                    Eigen::Matrix2d cov = Eigen::Matrix2d::Zero();
                    double mx = sx / denom, my = sy / denom;
                    for (int j = 0; j < np; ++j) {
                        double dx = data.x[ids[j]] - mx;
                        double dy = data.y[ids[j]] - my;
                        cov(0,0) += dx*dx; cov(0,1) += dx*dy;
                        cov(1,0) += dx*dy; cov(1,1) += dy*dy;
                    }
                    cov /= np;
                    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(cov);
                    auto ev = eig.eigenvalues();  // ascending order
                    stats(ci, ci_map["elongation"]) = (ev(0) > 1e-10)
                        ? (ev(1) / ev(0)) : std::numeric_limits<double>::quiet_NaN();
                }

                // max_cluster_frac
                if (has_clmol && np > 0) {
                    std::unordered_map<int,int> cl_cnt;
                    for (int mol : ids) cl_cnt[bm_data.cluster_per_molecule[mol]]++;
                    int mc = 0;
                    for (auto& [k, v] : cl_cnt) if (v > mc) mc = v;
                    stats(ci, ci_map["max_cluster_frac"]) = static_cast<double>(mc) / np;
                }

                // lifespan
                if (has_lspan) {
                    int guid = bm_data.components[ci].guid;
                    auto it2 = lifespan_map.find(guid);
                    stats(ci, ci_map["lifespan"]) = (it2 != lifespan_map.end()) ? it2->second : -1;
                }
            }

            save_cell_stat_df(stats, cell_names, col_names, out_paths.cell_stats);
        }

        // Save cell polygon hulls as GeoJSON
        if (polygon_format != "none") {
            spdlog::info("Saving cell polygons...");
            auto ids_by_cell_poly = split_ids(bm_data.assignment, n_cells_final, true);
            PolygonCollection polygons;
            polygons.reserve(n_cells_final);
            for (int ci = 0; ci < n_cells_final; ++ci) {
                const auto& ids = ids_by_cell_poly[ci];
                int np = static_cast<int>(ids.size());
                std::string cname = "cell_" + std::to_string(ci + 1);
                if (np >= 3) {
                    Eigen::MatrixXd pos2d(2, np);
                    for (int j = 0; j < np; ++j) {
                        pos2d(0, j) = data.x[ids[j]];
                        pos2d(1, j) = data.y[ids[j]];
                    }
                    polygons[cname] = convex_hull(pos2d);
                } else {
                    polygons[cname] = Eigen::MatrixXd(2, 0);
                }
            }
            save_polygons_geojson(polygons, out_paths.polygons_2d, polygon_format);
        }

        // Save count matrix
        spdlog::info("Saving count matrix...");
        {
            int ng = data.n_genes();
            std::vector<std::string> cell_names_cm(n_cells_final);
            for (int i = 0; i < n_cells_final; ++i)
                cell_names_cm[i] = "cell_" + std::to_string(i + 1);

            // Build sparse count matrix: n_cells x n_genes
            // Using composition_data which is 0-based internally
            std::vector<Eigen::Triplet<float>> trips;
            auto ids_by_cell_cm = split_ids(bm_data.assignment, n_cells_final, true);
            for (int ci = 0; ci < n_cells_final; ++ci) {
                std::unordered_map<int,float> cnt;
                for (int mol : ids_by_cell_cm[ci]) {
                    int g = bm_data.composition_data[mol];  // 0-based
                    if (g >= 0 && g < ng) cnt[g] += 1.0f;
                }
                for (auto& [g, v] : cnt) {
                    trips.emplace_back(ci, g, v);
                }
            }
            Eigen::SparseMatrix<float> cm(n_cells_final, ng);
            cm.setFromTriplets(trips.begin(), trips.end());

            if (count_matrix_format == "tsv") {
                Eigen::SparseMatrix<double> cm_d = cm.cast<double>();
                save_matrix_to_tsv(cm_d, data.gene_names, cell_names_cm, out_paths.counts);
            } else {
                save_matrix_to_loom(cm, data.gene_names, cell_names_cm, out_paths.counts);
            }
        }

        // Save parameters dump
        save_params_toml(opts, cli_cmd, out_paths.params_dump);

        spdlog::info("Results saved to '{}'", output);
    };

    if (data.is_3d()) {
        run_segmentation(std::integral_constant<int, 3>{});
    } else {
        run_segmentation(std::integral_constant<int, 2>{});
    }

    return 0;
}

// ============================================================================
// Subcommand: preview
// ============================================================================

int cmd_preview(
    const std::string& coordinates,
    RunOptions& opts,
    const std::string& output
) {
    spdlog::info("Loading data from '{}'...", coordinates);

    fill_and_check_data_options(opts.data);

    auto data = load_molecules(coordinates, opts.data);
    spdlog::info("Loaded {} transcripts, {} genes.", data.n_molecules(), data.n_genes());

    fill_and_check_plotting_options(opts.plotting, opts.data.min_molecules_per_cell, data.n_genes());

    // Confidence estimation — done once; reuse knn, adj_list, and noise_result
    // everywhere below instead of recomputing them.
    spdlog::info("Estimating noise level...");
    int nn_id = opts.data.confidence_nn_id;
    if (nn_id <= 0) nn_id = std::max(data.n_genes() / 10, 10);

    auto pos = data.position_matrix();
    auto knn = knn_parallel(pos, pos, nn_id + 1, true);

    std::vector<double> edge_lengths(data.n_molecules());
    for (int i = 0; i < data.n_molecules(); ++i) {
        int k = static_cast<int>(knn.distances[i].size());
        edge_lengths[i] = (k > nn_id) ? knn.distances[i][nn_id] : knn.distances[i].back();
    }

    auto adj_list    = build_molecule_graph(data, false);
    auto noise_result = fit_noise_probabilities(edge_lengths, adj_list, nullptr, 100, 0.005, true);

    data.confidence.resize(data.n_molecules());
    for (int i = 0; i < data.n_molecules(); ++i)
        data.confidence[i] = noise_result.assignment_probs(i, 0);

    spdlog::info("Done. Noise estimation complete.");

    // Gene composition colors
    spdlog::info("Estimating local colors...");
    int comp_k = opts.plotting.gene_composition_neighborhood;
    auto neighb_cm = neighborhood_count_matrix(pos, data.gene, comp_k, data.n_genes());

    // Log-transform (matching Julia: log(nzval * 10000 + 1e-5))
    for (int k = 0; k < neighb_cm.outerSize(); ++k) {
        for (Eigen::SparseMatrix<float>::InnerIterator it(neighb_cm, k); it; ++it) {
            it.valueRef() = static_cast<float>(std::log(it.value() * 10000.0f + 1e-5f));
        }
    }

    auto mol_vecs   = estimate_gene_vectors(neighb_cm, data.gene, 20, "ri", true);
    auto gene_colors = gene_composition_color_embedding(mol_vecs, data.confidence);
    spdlog::info("Done.");

    // Gene structure embedding (reuses adj_list already built above)
    spdlog::info("Estimating gene structure...");
    auto gene_structure = estimate_gene_structure_embedding(
        data.gene, data.gene_names, data.confidence, adj_list);
    spdlog::info("Done.");

    // Generate HTML report
    spdlog::info("Generating HTML report...");
    auto html = generate_preview_html(data, gene_colors, edge_lengths, noise_result, nn_id, &gene_structure);

    std::ofstream out_file(output);
    if (!out_file) {
        spdlog::error("Could not write to '{}'", output);
        return 1;
    }
    out_file << html;
    out_file.close();

    spdlog::info("Preview saved to '{}'", output);
    return 0;
}

// ============================================================================
// Subcommand: segfree
// ============================================================================

int cmd_segfree(
    const std::string& coordinates,
    RunOptions& opts,
    int k_neighbors,
    const std::string& output
) {
    spdlog::info("Loading data from '{}'...", coordinates);

    fill_and_check_data_options(opts.data);

    auto data = load_molecules(coordinates, opts.data);
    spdlog::info("Loaded {} transcripts, {} genes.", data.n_molecules(), data.n_genes());

    if (k_neighbors <= 0) {
        k_neighbors = default_param_value(
            "composition_neighborhood", opts.data.min_molecules_per_cell, -1, data.n_genes());
    }
    spdlog::info("Using k={} neighbors for NCV composition.", k_neighbors);

    // Neighborhood count matrix: n_genes × n_mols sparse
    spdlog::info("Estimating neighborhoods...");
    auto pos = data.position_matrix();
    auto neighb_cm = neighborhood_count_matrix(pos, data.gene, k_neighbors, data.n_genes());

    // Log-transform (matching Julia and preview pipeline)
    for (int k = 0; k < neighb_cm.outerSize(); ++k) {
        for (Eigen::SparseMatrix<float>::InnerIterator it(neighb_cm, k); it; ++it) {
            it.valueRef() = static_cast<float>(std::log(it.value() * 10000.0f + 1e-5f));
        }
    }

    // Per-molecule confidence (noise model)
    spdlog::info("Estimating molecule confidences...");
    append_confidence(data, opts.data.confidence_nn_id);

    // Gene vectors via randomized indexing, then UMAP color embedding
    spdlog::info("Estimating gene colors...");
    auto mol_vecs = estimate_gene_vectors(neighb_cm, data.gene, 20, "ri", true);
    auto gene_colors = gene_composition_color_embedding(mol_vecs, data.confidence);

    // Cell/NCV names: "V{i}" (1-based), matching Julia's get_cell_name(:ncv)
    int n = data.n_molecules();
    std::vector<std::string> ncv_names(n);
    for (int i = 0; i < n; ++i) ncv_names[i] = "V" + std::to_string(i + 1);

    // Save Loom file.
    // Julia stores neighb_cm transposed: n_mols × n_genes (rows = cells, cols = genes).
    spdlog::info("Saving results to '{}'...", output);
    Eigen::SparseMatrix<float> ncv_mat = neighb_cm.transpose();

    LoomColAttrs col_attrs;
    col_attrs["ncv_color"]  = gene_colors;  // vector<string>
    col_attrs["confidence"] = std::vector<double>(data.confidence.begin(),
                                                   data.confidence.end());

    save_matrix_to_loom(ncv_mat, data.gene_names, ncv_names, output, col_attrs);

    spdlog::info("Done.");
    return 0;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    CLI::App app{"Baysor — Bayesian cell segmentation of spatial transcriptomics data"};
    app.require_subcommand(1);

    // Pre-scan argv for -c/--config so we can load config before CLI11 registers
    // option defaults. This way config values serve as defaults and explicit CLI
    // flags override them (proper config-then-CLI precedence).
    std::string config_path;
    for (int i = 1; i < argc - 1; ++i) {
        std::string arg = argv[i];
        if (arg == "-c" || arg == "--config") {
            config_path = argv[i + 1];
            break;
        }
    }

    RunOptions opts;
    if (!config_path.empty()) {
        try {
            opts = load_config(config_path);
            spdlog::info("Loaded config from '{}'", config_path);
        } catch (const std::exception& e) {
            spdlog::error("Failed to load config '{}': {}", config_path, e.what());
            return 1;
        }
    }

    // ---- run ----
    auto* run = app.add_subcommand("run", "Run cell segmentation");

    std::string run_coordinates, run_prior_seg;
    std::string run_output = "segmentation.csv";
    std::string run_polygon_format = "FeatureCollection";
    std::string run_count_format = "loom";
    bool run_plot = false;
    bool run_skip_ncv_color = false;

    run->add_option("coordinates", run_coordinates,
        "CSV or Parquet file with coordinates of molecules and gene type")
        ->required();
    run->add_option("prior_segmentation", run_prior_seg,
        "Image/MAT file with segmentation mask, or ':column_name' for a column in the coordinates file");

    run->add_option("-c,--config", config_path,
        "TOML file with configuration");
    run->add_option("-x,--x-column", opts.data.x_col,
        "Name of x column (default: x)");
    run->add_option("-y,--y-column", opts.data.y_col,
        "Name of y column (default: y)");
    run->add_option("-z,--z-column", opts.data.z_col,
        "Name of z column (default: z)");
    run->add_option("-g,--gene-column", opts.data.gene_col,
        "Name of gene column (default: gene)");
    run->add_option("-m,--min-molecules-per-cell", opts.data.min_molecules_per_cell,
        "Minimal number of molecules for a cell to be considered as real");
    run->add_option("-s,--scale", opts.segmentation.scale,
        "Approximate cell radius. Sets estimate-scale-from-centers to false");
    run->add_option("--scale-std", opts.segmentation.scale_std,
        "Std of scale across cells. Number or 'N%' relative to scale (default: 25%)");
    run->add_option("--n-clusters", opts.segmentation.n_clusters,
        "Number of molecule clusters / major cell types (default: 4)");
    run->add_option("--prior-segmentation-confidence", opts.segmentation.prior_segmentation_confidence,
        "Confidence of prior segmentation results, in [0,1] (default: 0.2)");
    run->add_option("--min-molecules-per-gene", opts.data.min_molecules_per_gene,
        "Minimal number of molecules per gene (default: 1)");
    run->add_option("--exclude-genes", opts.data.exclude_genes,
        "Comma-separated list of genes or patterns to exclude (e.g. 'Blank*,MALAT1')");
    run->add_option("--nuclei-genes", opts.segmentation.nuclei_genes,
        "Comma-separated list of nuclei-specific genes");
    run->add_option("--cyto-genes", opts.segmentation.cyto_genes,
        "Comma-separated list of cytoplasm-specific genes");
    run->add_option("-o,--output", run_output,
        "Output file or directory (default: segmentation.csv)");
    run->add_option("--polygon-format", run_polygon_format,
        "Polygon output format: FeatureCollection, GeometryCollection, or none (default: FeatureCollection)");
    run->add_option("--count-matrix-format", run_count_format,
        "Count matrix format: loom or tsv (default: loom)");
    run->add_flag("-p,--plot", run_plot,
        "Save an HTML diagnostic plot");
    run->add_flag("--skip-ncv-color", run_skip_ncv_color,
        "Skip neighborhood composition color embedding to speed up development runs");
    run->add_flag("--force-2d", opts.data.force_2d,
        "Ignore z-column in the data");
    run->add_option("--iters", opts.segmentation.iters,
        "Maximum number of algorithm iterations (default: 500)");
    run->add_option("--tol", opts.segmentation.tol,
        "Convergence tolerance: stop when <tol fraction of molecules change assignment "
        "over 20 consecutive iterations. 0 = always run all --iters (default: 0.005)");
    run->add_option("--n-cells-init", opts.segmentation.n_cells_init,
        "Initial number of cells (default: auto)");
    run->add_option("--unassigned-prior-label", opts.segmentation.unassigned_prior_label,
        "Label for unassigned cells in prior segmentation (default: 0)");

    // ---- preview ----
    auto* preview = app.add_subcommand("preview", "Plot a dataset preview");

    std::string prev_coordinates;
    std::string prev_output = "preview.html";

    preview->add_option("coordinates", prev_coordinates,
        "CSV or Parquet file with coordinates of molecules and gene type")
        ->required();
    preview->add_option("-c,--config", config_path,
        "TOML file with configuration");
    preview->add_option("-x,--x-column", opts.data.x_col,
        "Name of x column (default: x)");
    preview->add_option("-y,--y-column", opts.data.y_col,
        "Name of y column (default: y)");
    preview->add_option("-z,--z-column", opts.data.z_col,
        "Name of z column (default: z)");
    preview->add_option("-g,--gene-column", opts.data.gene_col,
        "Name of gene column (default: gene)");
    preview->add_option("-m,--min-molecules-per-cell", opts.data.min_molecules_per_cell,
        "Minimal number of molecules for a cell to be considered as real");
    preview->add_option("-o,--output", prev_output,
        "Output file or directory (default: preview.html)");
    preview->add_flag("--force-2d", opts.data.force_2d,
        "Ignore z-column in the data");

    // ---- segfree ----
    auto* segfree = app.add_subcommand("segfree", "Extract Neighborhood Composition Vectors (NCVs)");

    std::string sf_coordinates;
    std::string sf_output = "ncvs.loom";
    int sf_k_neighbors = 0;

    segfree->add_option("coordinates", sf_coordinates,
        "CSV or Parquet file with coordinates of molecules and gene type")
        ->required();
    segfree->add_option("-c,--config", config_path,
        "TOML file with configuration");
    segfree->add_option("-x,--x-column", opts.data.x_col,
        "Name of x column (default: x)");
    segfree->add_option("-y,--y-column", opts.data.y_col,
        "Name of y column (default: y)");
    segfree->add_option("-z,--z-column", opts.data.z_col,
        "Name of z column (default: z)");
    segfree->add_option("-g,--gene-column", opts.data.gene_col,
        "Name of gene column (default: gene)");
    segfree->add_option("-m,--min-molecules-per-cell", opts.data.min_molecules_per_cell,
        "Minimal number of molecules for a cell to be considered as real");
    segfree->add_option("-k,--k-neighbors", sf_k_neighbors,
        "Number of neighbors for segmentation-free pseudo-cells (default: inferred)");
    segfree->add_option("-o,--output", sf_output,
        "Output file or directory (default: ncvs.loom)");
    segfree->add_flag("--force-2d", opts.data.force_2d,
        "Ignore z-column in the data");

    // ---- Parse ----
    CLI11_PARSE(app, argc, argv);

    // Reconstruct CLI command string for params dump
    std::string cli_cmd;
    for (int i = 0; i < argc; ++i) {
        if (i > 0) cli_cmd += " ";
        cli_cmd += argv[i];
    }

    // Dispatch
    if (run->parsed()) {
        if (opts.segmentation.scale > 0) {
            opts.segmentation.estimate_scale_from_centers = false;
        }
        if (run_prior_seg.empty() && opts.segmentation.scale <= 0) {
            spdlog::error("Either prior_segmentation or --scale must be provided.");
            return 1;
        }
        return cmd_run(run_coordinates, run_prior_seg, opts, run_output,
                       run_plot, run_skip_ncv_color, run_polygon_format, run_count_format, cli_cmd);
    }

    if (preview->parsed()) {
        return cmd_preview(prev_coordinates, opts, prev_output);
    }

    if (segfree->parsed()) {
        return cmd_segfree(sf_coordinates, opts, sf_k_neighbors, sf_output);
    }

    return 0;
}
