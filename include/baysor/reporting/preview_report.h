#pragma once

#include "baysor/data_loading/data.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/reporting/color_utils.h"
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace baysor {

/// Render all molecules as a PNG (base64-encoded), coloured by the given hex strings.
/// width_px controls the output width; height is derived from the data aspect ratio.
std::string render_scatter_png(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::string>& colors,
    int width_px = 6000
);

/// Render confidence scatter as a PNG (base64-encoded), coloured by a blue-orange gradient.
std::string render_confidence_png(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& confidence,
    int width_px = 6000
);

/// Generate Vega-Lite spec: noise estimation histogram + fitted PDFs
nlohmann::json vega_noise_histogram(
    const std::vector<double>& edge_lengths,
    const std::vector<double>& confidence,
    double signal_mu, double signal_sigma,
    double noise_mu, double noise_sigma,
    int nn_id = 0,
    int n_bins = 50
);

/// Generate Vega-Lite spec: gene frequency stacked bar chart
nlohmann::json vega_gene_frequency(
    const std::vector<int>& genes,
    const std::vector<double>& confidence,
    const std::vector<std::string>& gene_names
);

/// Generate Vega-Lite spec: gene structure scatter (UMAP of genes by spatial co-occurrence)
nlohmann::json vega_gene_structure(const GeneStructureEmbedding& emb);

/// Generate complete HTML preview report
std::string generate_preview_html(
    const MoleculeData& data,
    const std::vector<std::string>& gene_colors,
    const std::vector<double>& edge_lengths,
    const NoiseFitResult& noise_result,
    int confidence_nn_id,
    const GeneStructureEmbedding* gene_structure = nullptr
);

} // namespace baysor
