#include "baysor/segmentation/segmentation.h"

#include "baysor/data_loading/prior_segmentation.h"
#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/processing/bmm_algorithm/tracing.h"
#include "baysor/processing/data_processing/initialization.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/processing/utils/convex_hull.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/utils/general.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <utility>

#ifndef BAYSOR_VERSION
#define BAYSOR_VERSION "unknown"
#endif

#ifndef BAYSOR_BUILD_REVISION
#define BAYSOR_BUILD_REVISION "unknown"
#endif

namespace baysor {

namespace {

std::uint64_t splitmix64_finalizer(std::uint64_t value) noexcept {
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

unsigned int narrow_seed(std::uint64_t seed) noexcept {
    return static_cast<unsigned int>(seed) ^ static_cast<unsigned int>(seed >> 32U);
}

int narrow_signed_seed(std::uint64_t seed) noexcept {
    return static_cast<int>(narrow_seed(seed) & 0x7fffffffU);
}

class ScopedNativeThreads {
public:
    explicit ScopedNativeThreads(int requested_threads)
        : previous_max_threads_(omp_get_max_threads()) {
        if (requested_threads > 0) {
            omp_set_num_threads(requested_threads);
        }
        effective_threads_ = omp_get_max_threads();
    }

    ~ScopedNativeThreads() {
        omp_set_num_threads(previous_max_threads_);
    }

    ScopedNativeThreads(const ScopedNativeThreads&) = delete;
    ScopedNativeThreads& operator=(const ScopedNativeThreads&) = delete;

    [[nodiscard]] int effective_threads() const noexcept {
        return effective_threads_;
    }

private:
    int previous_max_threads_ = 1;
    int effective_threads_ = 1;
};

[[nodiscard]] bool cancellation_requested(const CancellationToken& cancellation) noexcept {
    return cancellation.is_cancellation_requested();
}

void validate_request(const SegmentationRequest& request) {
    if (request.molecules.path.empty()) {
        throw SegmentationError(
            SegmentationErrorCode::InvalidRequest,
            "The molecule input path must not be empty.");
    }
    if (request.execution.native_threads < 0) {
        throw SegmentationError(
            SegmentationErrorCode::InvalidRequest,
            "native_threads must be zero or positive.");
    }
    if (request.segmentation.iters < 0) {
        throw SegmentationError(
            SegmentationErrorCode::InvalidRequest,
            "The number of segmentation iterations must not be negative.");
    }
    if (!request.segmentation.nuclei_genes.empty() || !request.segmentation.cyto_genes.empty()) {
        throw SegmentationError(
            SegmentationErrorCode::InvalidRequest,
            "Nucleus/cytoplasm compartment segmentation is not implemented by the native operation.");
    }
}

int infer_initial_cell_count(
    const MoleculeData& data,
    int min_molecules_per_cell
) {
    int inferred = default_param_value(
        "n_cells_init", min_molecules_per_cell, data.n_molecules());

    if (data.prior_segmentation.empty()) return inferred;

    const int max_label = *std::max_element(
        data.prior_segmentation.begin(), data.prior_segmentation.end());
    if (max_label <= 0) return inferred;

    std::vector<int> segment_counts(max_label + 1, 0);
    int unassigned_molecules = 0;
    for (int label : data.prior_segmentation) {
        if (label > 0) {
            ++segment_counts[label];
        } else {
            ++unassigned_molecules;
        }
    }

    int active_prior_segments = 0;
    for (int label = 1; label <= max_label; ++label) {
        if (segment_counts[label] > 0) ++active_prior_segments;
    }
    if (active_prior_segments == 0) return inferred;

    constexpr double prior_segment_multiplier = 2.25;
    constexpr double unassigned_multiplier = 2.0;
    int prior_based = static_cast<int>(std::ceil(
        prior_segment_multiplier * static_cast<double>(active_prior_segments) +
        unassigned_multiplier * static_cast<double>(unassigned_molecules) /
            std::max(min_molecules_per_cell, 1)));
    prior_based = std::max(prior_based, active_prior_segments);
    return std::min(inferred, prior_based);
}

std::vector<std::string> make_cell_ids(int n_cells) {
    std::vector<std::string> cell_ids(n_cells);
    for (int i = 0; i < n_cells; ++i) {
        cell_ids[i] = "cell_" + std::to_string(i + 1);
    }
    return cell_ids;
}

template<int N>
CellStatistics build_cell_statistics(
    const MoleculeData& data,
    const BmmData<N>& bm_data,
    const std::vector<std::string>& cell_ids
) {
    const int n_cells = bm_data.n_components();
    const auto ids_by_cell = split_ids(bm_data.assignment, n_cells, true);

    std::unordered_map<int, int> lifespan_map;
    if (!bm_data.assignment_history.empty()) {
        lifespan_map = estimate_component_lifespan(bm_data.assignment_history);
    }

    const bool has_cluster = !bm_data.cluster_per_cell.empty();
    const bool has_assignment_confidence = !bm_data.assignment_confidence.empty();
    const bool has_molecule_clusters = !bm_data.cluster_per_molecule.empty();
    const bool has_lifespan = !lifespan_map.empty();

    CellStatistics result;
    result.cell_ids = cell_ids;
    result.columns = {"x", "y"};
    if constexpr (N == 3) result.columns.push_back("z");
    if (has_cluster) result.columns.push_back("cluster");
    result.columns.push_back("n_transcripts");
    result.columns.push_back("density");
    result.columns.push_back("elongation");
    result.columns.push_back("area");
    result.columns.push_back("avg_confidence");
    if (has_assignment_confidence) result.columns.push_back("avg_assignment_confidence");
    if (has_molecule_clusters) result.columns.push_back("max_cluster_frac");
    if (has_lifespan) result.columns.push_back("lifespan");

    std::unordered_map<std::string, int> column_index;
    for (int i = 0; i < static_cast<int>(result.columns.size()); ++i) {
        column_index[result.columns[i]] = i;
    }

    result.values.resize(n_cells, static_cast<int>(result.columns.size()));
    result.values.fill(std::numeric_limits<double>::quiet_NaN());

    for (int cell = 0; cell < n_cells; ++cell) {
        const auto& ids = ids_by_cell[cell];
        const int molecule_count = static_cast<int>(ids.size());
        const double denominator = molecule_count > 0 ? molecule_count : 1.0;

        double x_sum = 0.0;
        double y_sum = 0.0;
        double z_sum = 0.0;
        double confidence_sum = 0.0;
        double assignment_confidence_sum = 0.0;
        for (int molecule : ids) {
            x_sum += data.x[molecule];
            y_sum += data.y[molecule];
            if constexpr (N == 3) z_sum += data.z[molecule];
            confidence_sum += data.confidence[molecule];
            if (has_assignment_confidence) {
                assignment_confidence_sum += bm_data.assignment_confidence[molecule];
            }
        }

        result.values(cell, column_index.at("x")) = x_sum / denominator;
        result.values(cell, column_index.at("y")) = y_sum / denominator;
        if constexpr (N == 3) {
            result.values(cell, column_index.at("z")) = z_sum / denominator;
        }
        result.values(cell, column_index.at("n_transcripts")) = molecule_count;
        result.values(cell, column_index.at("avg_confidence")) = confidence_sum / denominator;
        if (has_assignment_confidence) {
            result.values(cell, column_index.at("avg_assignment_confidence")) =
                assignment_confidence_sum / denominator;
        }
        if (has_cluster) {
            result.values(cell, column_index.at("cluster")) = bm_data.cluster_per_cell[cell];
        }

        if (molecule_count > 2) {
            Eigen::MatrixXd positions(2, molecule_count);
            for (int i = 0; i < molecule_count; ++i) {
                positions(0, i) = data.x[ids[i]];
                positions(1, i) = data.y[ids[i]];
            }
            const auto hull = convex_hull(positions);
            const double area = polygon_area(hull);
            result.values(cell, column_index.at("area")) = area;
            result.values(cell, column_index.at("density")) = area > 0.0
                ? molecule_count / area
                : std::numeric_limits<double>::quiet_NaN();

            Eigen::Matrix2d covariance = Eigen::Matrix2d::Zero();
            const double mean_x = x_sum / denominator;
            const double mean_y = y_sum / denominator;
            for (int molecule : ids) {
                const double dx = data.x[molecule] - mean_x;
                const double dy = data.y[molecule] - mean_y;
                covariance(0, 0) += dx * dx;
                covariance(0, 1) += dx * dy;
                covariance(1, 0) += dx * dy;
                covariance(1, 1) += dy * dy;
            }
            covariance /= molecule_count;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigen(covariance);
            const auto eigenvalues = eigen.eigenvalues();
            result.values(cell, column_index.at("elongation")) = eigenvalues(0) > 1e-10
                ? eigenvalues(1) / eigenvalues(0)
                : std::numeric_limits<double>::quiet_NaN();
        }

        if (has_molecule_clusters && molecule_count > 0) {
            std::unordered_map<int, int> cluster_counts;
            for (int molecule : ids) ++cluster_counts[bm_data.cluster_per_molecule[molecule]];
            int maximum = 0;
            for (const auto& [cluster, count] : cluster_counts) {
                (void)cluster;
                maximum = std::max(maximum, count);
            }
            result.values(cell, column_index.at("max_cluster_frac")) =
                static_cast<double>(maximum) / molecule_count;
        }

        if (has_lifespan) {
            const auto lifespan = lifespan_map.find(bm_data.components[cell].guid);
            result.values(cell, column_index.at("lifespan")) =
                lifespan == lifespan_map.end() ? -1 : lifespan->second;
        }
    }

    return result;
}

template<int N>
CellByGeneCounts build_count_matrix(
    const MoleculeData& data,
    const BmmData<N>& bm_data,
    const std::vector<std::string>& cell_ids
) {
    CellByGeneCounts result;
    result.cell_ids = cell_ids;
    result.gene_names = data.gene_names;
    result.values.resize(bm_data.n_components(), data.n_genes());

    std::vector<Eigen::Triplet<float>> entries;
    const auto ids_by_cell = split_ids(
        bm_data.assignment, bm_data.n_components(), /*drop_zero=*/true);
    for (int cell = 0; cell < bm_data.n_components(); ++cell) {
        std::unordered_map<int, float> counts;
        for (int molecule : ids_by_cell[cell]) {
            const int gene = bm_data.composition_data[molecule];
            if (gene >= 0 && gene < data.n_genes()) counts[gene] += 1.0F;
        }
        for (const auto& [gene, count] : counts) {
            entries.emplace_back(cell, gene, count);
        }
    }
    result.values.setFromTriplets(entries.begin(), entries.end());
    return result;
}

ConfidenceEstimationDiagnostics make_confidence_diagnostics(
    const ConfidenceEstimationDetails& details
) {
    ConfidenceEstimationDiagnostics result;
    result.edge_lengths = details.edge_lengths;
    result.fit_differences = details.fit_result.diffs;
    result.neighbor_index = details.nn_id;
    result.signal_mean = details.fit_result.signal_mu;
    result.signal_standard_deviation = details.fit_result.signal_sigma;
    result.noise_mean = details.fit_result.noise_mu;
    result.noise_standard_deviation = details.fit_result.noise_sigma;
    return result;
}

NeighborhoodCompositionDiagnostics make_neighborhood_diagnostics(
    const NcvReportEmbedding& report
) {
    NeighborhoodCompositionDiagnostics result;
    result.sample_molecule_indices = report.sample_ids;
    result.sample_embedding_x = report.sample_umap_x;
    result.sample_embedding_y = report.sample_umap_y;
    result.chosen_confidence_threshold = report.chosen_threshold;
    result.anchor_count = report.anchor_count;
    return result;
}

template<int N>
SegmentationOutcome finish_segmentation(
    MoleculeData data,
    const AdjList& adjacency,
    const ResolvedSegmentationOptions& options,
    const SegmentationProducts& requested_products,
    const ConfidenceEstimationDetails& confidence_details,
    const std::optional<ClusteringResult>& clustering_result,
    const std::vector<int>& molecule_clusters,
    Xoshiro256pp& core_random_state,
    std::uint64_t core_seed,
    std::uint64_t neighborhood_seed,
    std::uint64_t diagnostics_seed,
    const NativeRunProvenance& provenance,
    const CancellationToken& cancellation
) {
    auto bm_data = initialize_bmm_data<N>(
        data,
        adjacency,
        options.segmentation.n_cells_init,
        options.segmentation.scale,
        options.segmentation.scale_std,
        options.segmentation.prior_segmentation_confidence,
        options.molecules.min_molecules_per_cell,
        /*verbose=*/false);
    if (!molecule_clusters.empty()) bm_data.cluster_per_molecule = molecule_clusters;

    if (cancellation_requested(cancellation)) return SegmentationCancelled{};

    const int history_depth = std::max(1, options.segmentation.iters / 10);
    bmm(
        bm_data,
        /*min_molecules_drop=*/2,
        options.segmentation.iters,
        history_depth,
        /*verbose=*/false,
        /*component_split_step=*/3,
        /*refine=*/true,
        /*freeze_composition=*/false,
        /*freeze_position=*/false,
        /*freeze_components=*/false,
        options.segmentation.tol,
        options.molecules.min_molecules_per_cell,
        &core_random_state,
        core_seed);

    if (cancellation_requested(cancellation)) return SegmentationCancelled{};

    SegmentationResult result;
    result.cell_ids = make_cell_ids(bm_data.n_components());
    result.cell_assignments = bm_data.assignment;
    result.assignment_confidence = bm_data.assignment_confidence;
    result.molecule_clusters = molecule_clusters;

    std::optional<NcvReportEmbedding> ncv_report;
    if (requested_products.neighborhood_composition_colors) {
        const auto positions = data.position_matrix();
        const int color_seed = narrow_signed_seed(neighborhood_seed);
        const auto ncv_model = fit_ncv_projected_model(
            positions,
            data.gene,
            data.n_genes(),
            data.confidence,
            options.neighborhood_composition.neighborhood_size,
            options.segmentation.cluster_basis_sample_size,
            /*n_components=*/20,
            /*include_full_projection=*/true,
            narrow_seed(neighborhood_seed));
        result.neighborhood_composition_colors =
            gene_composition_color_embedding_streaming(
                positions,
                data.gene,
                data.n_genes(),
                data.confidence,
                options.neighborhood_composition.neighborhood_size,
                options.segmentation.cluster_basis_sample_size,
                /*sample_size=*/20000,
                color_seed,
                /*n_pca_dims=*/10,
                options.segmentation.cluster_graph_k,
                &ncv_model);
        if (requested_products.diagnostics) {
            ncv_report = gene_composition_report_embedding_streaming(
                positions,
                data.gene,
                data.n_genes(),
                data.confidence,
                options.neighborhood_composition.neighborhood_size,
                options.segmentation.cluster_basis_sample_size,
                /*sample_size=*/20000,
                narrow_signed_seed(diagnostics_seed),
                /*n_pca_dims=*/10,
                options.segmentation.cluster_graph_k,
                &ncv_model);
        }
    }

    if (requested_products.cell_statistics) {
        result.cell_statistics = build_cell_statistics(data, bm_data, result.cell_ids);
    }
    if (requested_products.boundaries) {
        auto boundaries = boundary_polygons_auto(
            data.position_matrix(),
            bm_data.assignment,
            /*estimate_per_z=*/N == 3,
            &result.cell_ids,
            /*verbose=*/false,
            &core_random_state);
        result.boundaries_2d = std::move(boundaries.first);
        if constexpr (N == 3) result.boundaries_3d = std::move(boundaries.second);
    }
    if (requested_products.count_matrix) {
        result.count_matrix = build_count_matrix(data, bm_data, result.cell_ids);
    }

    if (requested_products.diagnostics) {
        SegmentationDiagnostics diagnostics;
        diagnostics.confidence_estimation = make_confidence_diagnostics(confidence_details);
        diagnostics.component_count_trace.reserve(bm_data.n_components_trace.size());
        for (const auto& snapshot : bm_data.n_components_trace) {
            diagnostics.component_count_trace.push_back({snapshot});
        }
        if (clustering_result) {
            diagnostics.molecule_clustering = MoleculeClusteringDiagnostics{
                clustering_result->diffs,
                clustering_result->change_fracs
            };
        }
        if (ncv_report) {
            diagnostics.neighborhood_composition =
                make_neighborhood_diagnostics(*ncv_report);
        }
        result.diagnostics = std::move(diagnostics);
    }

    result.resolved_options = options;
    result.provenance = provenance;
    result.produced_products.molecule_assignments = true;
    result.produced_products.molecule_confidence = !data.confidence.empty();
    result.produced_products.assignment_confidence = !result.assignment_confidence.empty();
    result.produced_products.molecule_clusters = !result.molecule_clusters.empty();
    result.produced_products.neighborhood_composition_colors =
        !result.neighborhood_composition_colors.empty();
    result.produced_products.cell_statistics = result.cell_statistics.has_value();
    result.produced_products.boundaries = result.boundaries_2d.has_value();
    result.produced_products.count_matrix = result.count_matrix.has_value();
    result.produced_products.diagnostics = result.diagnostics.has_value();
    result.molecules = std::move(data);

    if (cancellation_requested(cancellation)) return SegmentationCancelled{};
    return result;
}

SegmentationOutcome run_segmentation_impl(
    const SegmentationRequest& request,
    const CancellationToken& cancellation
) {
    validate_request(request);
    if (cancellation_requested(cancellation)) return SegmentationCancelled{};

    ScopedNativeThreads thread_scope(request.execution.native_threads);

    ResolvedSegmentationOptions options;
    options.molecules = request.molecules.options;
    options.prior = request.prior;
    options.segmentation = request.segmentation;
    options.neighborhood_composition = request.neighborhood_composition;
    options.execution.native_threads = thread_scope.effective_threads();

    try {
        fill_and_check_molecule_input_options(options.molecules);
        if (options.segmentation.scale > 0.0) {
            options.prior.estimate_scale_from_prior = false;
        }
        fill_and_check_prior_input_options(
            options.prior, options.molecules.min_molecules_per_cell);
        if (options.segmentation.n_clusters <= 0) {
            options.segmentation.n_clusters =
                default_cluster_count(options.segmentation.cluster_method);
        }
    } catch (const std::exception& error) {
        throw SegmentationError(SegmentationErrorCode::InvalidRequest, error.what());
    }

    if (options.prior.type == PriorInputType::None && options.segmentation.scale <= 0.0) {
        throw SegmentationError(
            SegmentationErrorCode::InvalidRequest,
            "Scale cannot be determined without a prior segmentation; provide a positive scale.");
    }

    MoleculeData data;
    try {
        data = load_molecules(request.molecules.path, options.molecules, options.prior);
    } catch (const std::exception& error) {
        throw SegmentationError(SegmentationErrorCode::MoleculeInput, error.what());
    }
    if (data.n_molecules() == 0) {
        throw SegmentationError(
            SegmentationErrorCode::MoleculeInput,
            "No molecules remain after loading and filtering the input.");
    }
    if (cancellation_requested(cancellation)) return SegmentationCancelled{};

    try {
        PlottingOptions plotting;
        plotting.gene_composition_neighborhood =
            options.neighborhood_composition.neighborhood_size;
        plotting.ncv_method = options.neighborhood_composition.method;
        fill_and_check_plotting_options(
            plotting, options.molecules.min_molecules_per_cell, data.n_genes());
        options.neighborhood_composition.neighborhood_size =
            plotting.gene_composition_neighborhood;
        options.neighborhood_composition.method = plotting.ncv_method;
    } catch (const std::exception& error) {
        throw SegmentationError(SegmentationErrorCode::InvalidRequest, error.what());
    }

    if (options.prior.type != PriorInputType::None) {
        try {
            const auto [scale, scale_standard_deviation] = load_prior_segmentation(
                data, options.prior, options.molecules.min_molecules_per_cell);
            if (options.prior.estimate_scale_from_prior && scale > 0.0) {
                options.segmentation.scale = scale;
                options.segmentation.scale_std = std::to_string(scale_standard_deviation);
            }
        } catch (const std::exception& error) {
            throw SegmentationError(SegmentationErrorCode::PriorInput, error.what());
        }
    }

    if (options.segmentation.n_cells_init <= 0) {
        options.segmentation.n_cells_init = infer_initial_cell_count(
            data, options.molecules.min_molecules_per_cell);
    }
    if (options.segmentation.scale <= 0.0) {
        throw SegmentationError(
            SegmentationErrorCode::PriorInput,
            "The segmentation scale could not be determined from the supplied prior.");
    }

    const auto core_seed = derive_random_substream_seed(
        request.random_seed, RandomSubstream::CoreScientific);
    const auto clustering_seed = derive_random_substream_seed(
        request.random_seed, RandomSubstream::MoleculeClustering);
    const auto neighborhood_seed = derive_random_substream_seed(
        request.random_seed, RandomSubstream::NeighborhoodComposition);
    const auto diagnostics_seed = derive_random_substream_seed(
        request.random_seed, RandomSubstream::Diagnostics);
    Xoshiro256pp core_random_state(core_seed);

    const auto confidence_details = estimate_confidence_details(
        data,
        options.molecules.confidence_nn_id,
        options.segmentation.prior_segmentation_confidence,
        &core_random_state,
        /*verbose=*/false);
    data.confidence.resize(data.n_molecules());
    for (int i = 0; i < data.n_molecules(); ++i) {
        data.confidence[i] = confidence_details.fit_result.assignment_probs(i, 0);
    }

    const auto adjacency = build_molecule_graph(
        data,
        /*filter=*/true,
        /*use_local_gene_similarities=*/false,
        AdjacencyType::Auto,
        /*composition_neighborhood=*/0,
        /*n_gene_pcs=*/0,
        &core_random_state);
    if (cancellation_requested(cancellation)) return SegmentationCancelled{};

    ClusteringOptions clustering_options;
    clustering_options.method = options.segmentation.cluster_method;
    clustering_options.n_clusters = options.segmentation.n_clusters;
    clustering_options.resolution = options.segmentation.cluster_resolution;
    clustering_options.graph_k = options.segmentation.cluster_graph_k;
    clustering_options.spatial_k = options.neighborhood_composition.neighborhood_size;
    clustering_options.n_dims = options.segmentation.cluster_n_dims;
    clustering_options.basis_sample_size = options.segmentation.cluster_basis_sample_size;
    clustering_options.random_seed = narrow_seed(clustering_seed);

    const bool run_clustering =
        clustering_options.method == ClusterMethod::Louvain ||
        clustering_options.method == ClusterMethod::Leiden ||
        (clustering_options.method == ClusterMethod::Mrf && clustering_options.n_clusters > 1);
    std::optional<ClusteringResult> clustering_result;
    std::vector<int> molecule_clusters;
    if (run_clustering) {
        clustering_result = cluster_molecules(
            clustering_options.method == ClusterMethod::Louvain ||
                    clustering_options.method == ClusterMethod::Leiden
                ? data.position_matrix()
                : Eigen::MatrixXd(),
            data.gene,
            adjacency,
            data.confidence,
            clustering_options,
            /*verbose=*/false);
        molecule_clusters = clustering_result->assignment;
    }

    NativeRunProvenance provenance;
    provenance.baysor_version = BAYSOR_VERSION;
    provenance.build_revision = BAYSOR_BUILD_REVISION;
    provenance.random_seed = request.random_seed;
    provenance.random_substreams = {
        {RandomSubstream::CoreScientific, core_seed},
        {RandomSubstream::MoleculeClustering, clustering_seed},
        {RandomSubstream::NeighborhoodComposition, neighborhood_seed},
        {RandomSubstream::Diagnostics, diagnostics_seed}
    };
    provenance.effective_native_threads = thread_scope.effective_threads();

    if (data.is_3d()) {
        return finish_segmentation<3>(
            std::move(data), adjacency, options, request.requested_products,
            confidence_details, clustering_result, molecule_clusters,
            core_random_state, core_seed, neighborhood_seed, diagnostics_seed,
            provenance, cancellation);
    }
    return finish_segmentation<2>(
        std::move(data), adjacency, options, request.requested_products,
        confidence_details, clustering_result, molecule_clusters,
        core_random_state, core_seed, neighborhood_seed, diagnostics_seed,
        provenance, cancellation);
}

} // namespace

std::uint64_t derive_random_substream_seed(
    std::uint64_t master_seed,
    RandomSubstream stream
) noexcept {
    if (master_seed == kDefaultSegmentationSeed) {
        return stream == RandomSubstream::CoreScientific ? 1ULL : 42ULL;
    }

    constexpr std::uint64_t golden_ratio = 0x9e3779b97f4a7c15ULL;
    const auto stream_id = static_cast<std::uint64_t>(stream);
    return splitmix64_finalizer(master_seed + golden_ratio * (stream_id + 1ULL));
}

const char* random_substream_name(RandomSubstream stream) noexcept {
    switch (stream) {
        case RandomSubstream::CoreScientific:
            return "core_scientific";
        case RandomSubstream::MoleculeClustering:
            return "molecule_clustering";
        case RandomSubstream::NeighborhoodComposition:
            return "neighborhood_composition";
        case RandomSubstream::Diagnostics:
            return "diagnostics";
    }
    return "unknown";
}

SegmentationError::SegmentationError(SegmentationErrorCode code, std::string message)
    : std::runtime_error(std::move(message)), code_(code) {}

SegmentationErrorCode SegmentationError::code() const noexcept {
    return code_;
}

SegmentationOutcome run_segmentation(
    const SegmentationRequest& request,
    const CancellationToken& cancellation
) {
    try {
        return run_segmentation_impl(request, cancellation);
    } catch (const SegmentationError&) {
        throw;
    } catch (const std::exception& error) {
        throw SegmentationError(SegmentationErrorCode::NativeProcessing, error.what());
    }
}

} // namespace baysor
