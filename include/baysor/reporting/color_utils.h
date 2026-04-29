#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace baysor {

struct AdjList;

struct NcvBasisModel {
    std::vector<int> basis_ids;
    int spatial_k = 0;
    double distance_floor = -1.0;
    Eigen::MatrixXf gene_emb_t;
    Eigen::MatrixXf basis_vecs;
};

struct NcvProjectedModel {
    NcvBasisModel basis;
    Eigen::MatrixXf mol_vecs;
};

/// Normalize a 3 x N embedding to LAB color range.
/// Matches Julia's normalize_embedding_to_lab_range!
void normalize_embedding_to_lab_range(
    Eigen::MatrixXd& embedding,
    double l_min = 10.0, double l_max = 90.0,
    double trim_frac = 0.0125,
    bool log_colors = false
);

/// Convert a 3 x N LAB embedding to hex color strings.
/// Row 0 = L, Row 1 = a, Row 2 = b.
std::vector<std::string> embedding_to_hex(const Eigen::MatrixXd& lab_embedding);

/// UMAP-based color embedding: fit 3-component UMAP on high-confidence sample,
/// interpolate all molecules, normalize to LAB, return hex colors.
/// mol_vecs: n_components x n_molecules matrix (float; kept sparse-friendly).
/// confidence: per-molecule confidence (for selecting the sample).
/// n_pca_dims: number of PCA dimensions used for the KNN lookup during
///   interpolation. Higher values preserve more cell-type structure at the
///   cost of slightly slower KD-tree queries (sweet spot: 10–15).
std::vector<std::string> gene_composition_color_embedding(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int sample_size = 20000,
    int seed = 42,
    int n_pca_dims = 10,
    int graph_k = 15
);

/// Report-oriented NCV embedding diagnostics.
/// sample_ids indexes the anchor molecules used to fit the NCV embedding.
/// sample_umap_* stores a 2D UMAP of those same anchors for report visualization.
struct NcvReportEmbedding {
    std::vector<std::string> colors;
    std::vector<int> sample_ids;
    std::vector<double> sample_umap_x;
    std::vector<double> sample_umap_y;
    double chosen_threshold = 0.95;
    int anchor_count = 0;
};

/// UMAP-based color embedding with extra diagnostics for `run --plot`.
/// This computes the same NCV colors as gene_composition_color_embedding, plus
/// a 2D UMAP on the sampled anchor molecules for report visualization.
NcvReportEmbedding gene_composition_report_embedding(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int sample_size = 20000,
    int seed = 42,
    int n_pca_dims = 10,
    int graph_k = 15
);

NcvBasisModel fit_ncv_basis_model(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int n_genes,
    const std::vector<double>& confidence,
    int k_neighbors,
    int basis_sample_size = 100000,
    int n_components = 20
);

NcvProjectedModel fit_ncv_projected_model(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int n_genes,
    const std::vector<double>& confidence,
    int k_neighbors,
    int basis_sample_size = 100000,
    int n_components = 20,
    bool include_full_projection = true
);

std::vector<std::string> gene_composition_color_embedding_streaming(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int n_genes,
    const std::vector<double>& confidence,
    int k_neighbors,
    int basis_sample_size = 100000,
    int sample_size = 20000,
    int seed = 42,
    int n_pca_dims = 10,
    int graph_k = 15,
    const NcvProjectedModel* precomputed_model = nullptr
);

NcvReportEmbedding gene_composition_report_embedding_streaming(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int n_genes,
    const std::vector<double>& confidence,
    int k_neighbors,
    int basis_sample_size = 100000,
    int sample_size = 20000,
    int seed = 42,
    int n_pca_dims = 10,
    int graph_k = 15,
    const NcvProjectedModel* precomputed_model = nullptr
);

/// Pairwise gene spatial co-occurrence matrix from the molecule adjacency graph.
/// Returns n_genes x n_genes matrix normalized by sqrt of sum weights.
Eigen::MatrixXd pairwise_gene_spatial_cor(
    const std::vector<int>& genes,
    const std::vector<double>& confidence,
    const AdjList& adj_list,
    double confidence_threshold = 0.95
);

/// 2D UMAP embedding of genes based on spatial co-occurrence.
struct GeneStructureEmbedding {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<std::string> gene_names;
    std::vector<double> marker_sizes;  // log(molecule_count) per gene
};

GeneStructureEmbedding estimate_gene_structure_embedding(
    const std::vector<int>& genes,
    const std::vector<std::string>& gene_names,
    const std::vector<double>& confidence,
    const AdjList& adj_list,
    int seed = 42
);

} // namespace baysor
