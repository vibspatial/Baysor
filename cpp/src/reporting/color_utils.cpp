#include "baysor/reporting/color_utils.h"
#include "baysor/processing/data_processing/umap_wrappers.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/utils/utils.h"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <cstdio>
#include <omp.h>

namespace baysor {

// ============================================================================
// Quantile helper (for internal use)
// ============================================================================

static double quantile_vec(const std::vector<double>& v, double p) {
    if (v.empty()) return 0.0;
    std::vector<double> s = v;
    std::sort(s.begin(), s.end());
    double idx = p * (s.size() - 1);
    int lo = static_cast<int>(std::floor(idx));
    int hi = std::min(lo + 1, static_cast<int>(s.size()) - 1);
    double frac = idx - lo;
    return s[lo] * (1.0 - frac) + s[hi] * frac;
}

// ============================================================================
// normalize_embedding_to_lab_range
// ============================================================================

void normalize_embedding_to_lab_range(
    Eigen::MatrixXd& embedding,
    double l_min, double l_max,
    double trim_frac,
    bool log_colors
) {
    int n = static_cast<int>(embedding.cols());
    if (n == 0) return;

    // Per-row: subtract quantile(row, trim_frac), clamp to 0
    for (int r = 0; r < 3; ++r) {
        std::vector<double> row_vals(n);
        for (int i = 0; i < n; ++i) row_vals[i] = embedding(r, i);
        double q_lo = quantile_vec(row_vals, trim_frac);
        for (int i = 0; i < n; ++i) {
            embedding(r, i) = std::max(embedding(r, i) - q_lo, 0.0);
        }
    }

    // Global max quantile
    std::vector<double> all_vals;
    all_vals.reserve(3 * n);
    for (int r = 0; r < 3; ++r) {
        for (int i = 0; i < n; ++i) {
            all_vals.push_back(embedding(r, i));
        }
    }
    double max_val = quantile_vec(all_vals, 1.0 - trim_frac);
    if (max_val > 0) {
        embedding /= max_val;
    }
    embedding = embedding.cwiseMin(1.0);

    if (log_colors) {
        // Compute global 5th percentile
        all_vals.clear();
        for (int r = 0; r < 3; ++r) {
            for (int i = 0; i < n; ++i) {
                all_vals.push_back(embedding(r, i));
            }
        }
        double q05 = std::max(quantile_vec(all_vals, 0.05), 1e-3);
        for (int r = 0; r < 3; ++r) {
            for (int i = 0; i < n; ++i) {
                embedding(r, i) = std::log10(embedding(r, i) + q05);
            }
        }
        // Per-row normalize to [0, 1]
        for (int r = 0; r < 3; ++r) {
            double mn = embedding.row(r).minCoeff();
            embedding.row(r).array() -= mn;
            double mx = embedding.row(r).maxCoeff();
            if (mx > 0) embedding.row(r) /= mx;
        }
    }

    // L channel: scale to [l_min, l_max]
    embedding.row(0) *= (l_max - l_min);
    embedding.row(0).array() += l_min;

    // a, b channels: scale to [-100, 100]
    embedding.row(1).array() -= 0.5;
    embedding.row(1) *= 200.0;
    embedding.row(2).array() -= 0.5;
    embedding.row(2) *= 200.0;
}

// ============================================================================
// LAB → sRGB → hex conversion
// ============================================================================

// Standard LAB → XYZ → linear RGB → sRGB pipeline
// Reference white: D65 (X=0.95047, Y=1.0, Z=1.08883)

static inline double lab_f_inv(double t) {
    constexpr double delta = 6.0 / 29.0;
    if (t > delta) return t * t * t;
    return 3.0 * delta * delta * (t - 4.0 / 29.0);
}

static inline double srgb_gamma(double c) {
    if (c <= 0.0031308) return 12.92 * c;
    return 1.055 * std::pow(c, 1.0 / 2.4) - 0.055;
}

static inline int clamp_byte(double v) {
    int b = static_cast<int>(std::round(v * 255.0));
    return std::max(0, std::min(255, b));
}

static std::string lab_to_hex(double L, double a, double b) {
    // LAB → XYZ
    double fy = (L + 16.0) / 116.0;
    double fx = a / 500.0 + fy;
    double fz = fy - b / 200.0;

    double X = 0.95047 * lab_f_inv(fx);
    double Y = 1.00000 * lab_f_inv(fy);
    double Z = 1.08883 * lab_f_inv(fz);

    // XYZ → linear RGB (sRGB D65)
    double rlin =  3.2404542 * X - 1.5371385 * Y - 0.4985314 * Z;
    double glin = -0.9692660 * X + 1.8760108 * Y + 0.0415560 * Z;
    double blin =  0.0556434 * X - 0.2040259 * Y + 1.0572252 * Z;

    // Linear RGB → sRGB with gamma
    int R = clamp_byte(srgb_gamma(rlin));
    int G = clamp_byte(srgb_gamma(glin));
    int B = clamp_byte(srgb_gamma(blin));

    char buf[8];
    std::snprintf(buf, sizeof(buf), "#%02X%02X%02X", R, G, B);
    return std::string(buf);
}

std::vector<std::string> embedding_to_hex(const Eigen::MatrixXd& lab_embedding) {
    int n = static_cast<int>(lab_embedding.cols());
    std::vector<std::string> colors(n);
    for (int i = 0; i < n; ++i) {
        colors[i] = lab_to_hex(lab_embedding(0, i), lab_embedding(1, i), lab_embedding(2, i));
    }
    return colors;
}

// ============================================================================
// gene_composition_color_embedding (UMAP-based)
// ============================================================================

std::vector<std::string> gene_composition_color_embedding(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int sample_size,
    int seed,
    int n_pca_dims
) {
    int n_components = static_cast<int>(mol_vecs.rows());
    int n_mols       = static_cast<int>(mol_vecs.cols());

    if (n_mols == 0 || n_components < 3) {
        return std::vector<std::string>(n_mols, "#808080");
    }

    // Select a spatially-uniform high-confidence sample for fitting the UMAP.
    // Matches Julia's select_ids_uniformly: sum ALL component dimensions (not just 0+1)
    // to get a spread measure, then pick evenly-spaced indices from the sorted order.
    sample_size = std::min(sample_size, n_mols);

    std::vector<int> high_conf_ids;
    high_conf_ids.reserve(n_mols);
    for (int i = 0; i < n_mols; ++i) {
        if (confidence[i] >= 0.95) high_conf_ids.push_back(i);
    }
    // Julia fallback: use only high_conf_ids (not all molecules) when count < sample_size
    if (static_cast<int>(high_conf_ids.size()) < sample_size) {
        sample_size = static_cast<int>(high_conf_ids.size());
    }

    // Sort by sum of all component dimensions — matches Julia's sum(vals, dims=2).
    // Precompute sums to avoid recomputing them O(N log N) times in the comparator.
    int hc = static_cast<int>(high_conf_ids.size());
    std::vector<std::pair<float, int>> sum_ids(hc);
    for (int i = 0; i < hc; ++i) {
        float s = mol_vecs.col(high_conf_ids[i]).sum();
        sum_ids[i] = {s, high_conf_ids[i]};
    }
    std::sort(sum_ids.begin(), sum_ids.end());
    for (int i = 0; i < hc; ++i) high_conf_ids[i] = sum_ids[i].second;

    std::vector<int> sample_ids;
    sample_ids.reserve(sample_size);
    for (int i = 0; i < sample_size; ++i) {
        int idx = static_cast<int>(std::round(static_cast<double>(i) * (hc - 1) / (sample_size - 1)));
        sample_ids.push_back(high_conf_ids[idx]);
    }

    // Build float sample matrix (n_components x sample_size) — no cast needed.
    Eigen::MatrixXf sample_mat(n_components, sample_size);
    for (int i = 0; i < sample_size; ++i)
        sample_mat.col(i) = mol_vecs.col(sample_ids[i]);

    // UMAP fit: use spread=2.0 to match Julia's UmapFit defaults, which produce
    // better colour separation than the umappp default of spread=1.0.
    Eigen::MatrixXd sample_emb = umap_embed(sample_mat.cast<double>(), 3,
        /*n_neighbors=*/15, /*n_epochs=*/200, seed, /*spread=*/2.0);

    // PCA-reduce sample from n_components → n_pca dimensions before building the
    // interpolation KNN index.  In 20D the VPtree degenerates to exhaustive search
    // (curse of dimensionality), making interpolation O(N_mols × N_sample).
    // In 3D the VPtree search is O(log N_sample), giving a ~1000× speedup.
    // The first few PCs capture the dominant variance in gene-expression space,
    // so neighbour quality changes are minimal for colour interpolation.
    const int n_pca = std::min(n_pca_dims, n_components);
    Eigen::VectorXf sample_mean = sample_mat.rowwise().mean();
    Eigen::MatrixXf sample_centered = sample_mat.colwise() - sample_mean;

    // Thin U of (n_components × sample_size): columns are principal components.
    Eigen::BDCSVD<Eigen::MatrixXf> svd(sample_centered, Eigen::ComputeThinU);
    Eigen::MatrixXf pca_basis = svd.matrixU().leftCols(n_pca);  // n_components × n_pca

    // Project sample to PCA space: n_pca × sample_size
    Eigen::MatrixXf sample_pca =
        (pca_basis.transpose() * sample_centered).eval();

    // Project all molecules to PCA space: n_pca × n_mols
    Eigen::MatrixXd all_pca =
        (pca_basis.transpose() * (mol_vecs.colwise() - sample_mean))
        .cast<double>();

    // KNN in PCA-3D space: nanoflann KD-tree via knn_parallel (already OMP-parallel).
    // tree = sample_pca (n_pca × sample_size), query = all_pca (n_pca × n_mols).
    int k_interp = std::min(5, sample_size - 1);
    auto knn = knn_parallel(sample_pca.cast<double>(), all_pca, k_interp);

    // Weighted interpolation of UMAP coordinates.
    Eigen::MatrixXd emb(3, n_mols);
    constexpr double dist_offset = 1e-10;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_mols; ++i) {
        double w_sum = 0.0;
        Eigen::Vector3d weighted = Eigen::Vector3d::Zero();
        for (int j = 0; j < static_cast<int>(knn.indices[i].size()); ++j) {
            double w = 1.0 / (knn.distances[i][j] + dist_offset);
            weighted += w * sample_emb.col(knn.indices[i][j]);
            w_sum    += w;
        }
        emb.col(i) = weighted / w_sum;
    }

    normalize_embedding_to_lab_range(emb);
    return embedding_to_hex(emb);
}

// ============================================================================
// pairwise_gene_spatial_cor
// ============================================================================

Eigen::MatrixXd pairwise_gene_spatial_cor(
    const std::vector<int>& genes,
    const std::vector<double>& confidence,
    const AdjList& adj_list,
    double confidence_threshold
) {
    int n = static_cast<int>(genes.size());
    int n_genes = *std::max_element(genes.begin(), genes.end());  // 1-based max

    Eigen::MatrixXd cor_mat = Eigen::MatrixXd::Zero(n_genes, n_genes);
    std::vector<double> sum_weight(n_genes, 0.0);

    // Accumulate co-occurrence counts in per-thread local buffers to avoid
    // false-sharing and atomic overhead, then reduce at the end.
    #pragma omp parallel
    {
        Eigen::MatrixXd local_cor = Eigen::MatrixXd::Zero(n_genes, n_genes);
        std::vector<double> local_sw(n_genes, 0.0);

        #pragma omp for schedule(dynamic, 1024)
        for (int gi = 0; gi < n; ++gi) {
            if (confidence[gi] < confidence_threshold) continue;
            int g2 = genes[gi] - 1;  // 0-based
            int nc = adj_list.neighbor_count(gi);
            const int32_t* nb_ids = adj_list.neighbor_ids(gi);
            const double*  nb_wts = adj_list.neighbor_weights(gi);

            for (int ai = 0; ai < nc; ++ai) {
                int nb = nb_ids[ai];
                if (confidence[nb] < confidence_threshold) continue;
                int g1 = genes[nb] - 1;  // 0-based
                double cw = nb_wts[ai];
                local_cor(g2, g1) += cw;
                local_sw[g1]      += cw;
                local_sw[g2]      += cw;
            }
        }

        #pragma omp critical
        {
            cor_mat    += local_cor;
            for (int i = 0; i < n_genes; ++i) sum_weight[i] += local_sw[i];
        }
    }

    for (int ci = 0; ci < n_genes; ++ci) {
        for (int ri = 0; ri < n_genes; ++ri) {
            double denom = std::sqrt(sum_weight[ri] * sum_weight[ci]);
            cor_mat(ri, ci) /= std::max(denom, 0.1);
        }
    }

    return cor_mat;
}

// ============================================================================
// estimate_gene_structure_embedding
// ============================================================================

GeneStructureEmbedding estimate_gene_structure_embedding(
    const std::vector<int>& genes,
    const std::vector<std::string>& gene_names,
    const std::vector<double>& confidence,
    const AdjList& adj_list,
    int seed
) {
    int n_genes = static_cast<int>(gene_names.size());

    Eigen::MatrixXd cor_mat = pairwise_gene_spatial_cor(genes, confidence, adj_list);

    // Convert correlation to distance, matching Julia:
    //   p_dists = 1 - clamp(cor, min_cor, max_cor) / max_cor
    std::vector<double> flat(cor_mat.data(), cor_mat.data() + cor_mat.size());
    std::vector<double> pos_vals;
    for (double v : flat) if (v > 0) pos_vals.push_back(v);
    std::sort(pos_vals.begin(), pos_vals.end());

    double min_cor = 0.0, max_cor = 1.0;
    if (!pos_vals.empty()) {
        int idx001 = static_cast<int>(0.001 * (pos_vals.size() - 1));
        int idx99  = static_cast<int>(0.99  * (flat.size()    - 1));
        std::nth_element(flat.begin(), flat.begin() + idx99, flat.end());
        min_cor = pos_vals[idx001] / 5.0;
        max_cor = std::max(flat[idx99], 1e-10);
    }

    Eigen::MatrixXd dist_mat(n_genes, n_genes);
    for (int ci = 0; ci < n_genes; ++ci) {
        for (int ri = 0; ri < n_genes; ++ri) {
            double c = std::max(std::min(cor_mat(ri, ci), max_cor), min_cor);
            dist_mat(ri, ci) = 1.0 - c / max_cor;
        }
    }
    dist_mat.diagonal().setZero();

    int k = std::max(std::min(15, n_genes / 2), 2);
    Eigen::MatrixXd emb = umap_embed_precomputed(dist_mat, 2, k, 5000, seed);

    // Molecule counts per gene (for marker size, matching Julia's log(count_array)).
    std::vector<int> counts(n_genes, 0);
    for (int g : genes) {
        int gi = g - 1;
        if (gi >= 0 && gi < n_genes) counts[gi]++;
    }

    GeneStructureEmbedding result;
    result.x.resize(n_genes);
    result.y.resize(n_genes);
    result.gene_names = gene_names;
    result.marker_sizes.resize(n_genes);
    for (int i = 0; i < n_genes; ++i) {
        result.x[i]            = emb(0, i);
        result.y[i]            = emb(1, i);
        result.marker_sizes[i] = std::log(static_cast<double>(std::max(counts[i], 1)));
    }
    return result;
}

} // namespace baysor
