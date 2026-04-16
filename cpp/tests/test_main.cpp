#include <gtest/gtest.h>

#include "baysor/data_loading/data.h"
#include "baysor/data_loading/prior_segmentation.h"
#include "baysor/utils/general.h"
#include "baysor/utils/options.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/processing/data_processing/triangulation.h"
#include "baysor/processing/data_processing/initialization.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/reporting/color_utils.h"

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>

// Helper: write a temp CSV file and return its path
static std::string write_temp_csv(const std::string& content, const std::string& suffix = ".csv") {
    char tmpl[] = "/tmp/baysor_test_XXXXXX";
    int fd = mkstemp(tmpl);
    EXPECT_GE(fd, 0);
    close(fd);
    std::string path = std::string(tmpl) + suffix;
    std::rename(tmpl, path.c_str());
    std::ofstream f(path);
    f << content;
    f.close();
    return path;
}

// ============================================================================
// Gene encoding
// ============================================================================

TEST(GeneEncoding, BasicEncoding) {
    baysor::MoleculeData data;
    std::vector<std::string> genes = {"GeneC", "GeneA", "GeneB", "GeneA", "GeneC"};
    baysor::encode_genes(data, genes);

    // gene_names should be sorted alphabetically
    ASSERT_EQ(data.gene_names.size(), 3u);
    EXPECT_EQ(data.gene_names[0], "GeneA");
    EXPECT_EQ(data.gene_names[1], "GeneB");
    EXPECT_EQ(data.gene_names[2], "GeneC");

    // gene IDs should be 1-based
    ASSERT_EQ(data.gene.size(), 5u);
    EXPECT_EQ(data.gene[0], 3);  // GeneC -> 3
    EXPECT_EQ(data.gene[1], 1);  // GeneA -> 1
    EXPECT_EQ(data.gene[2], 2);  // GeneB -> 2
    EXPECT_EQ(data.gene[3], 1);  // GeneA -> 1
    EXPECT_EQ(data.gene[4], 3);  // GeneC -> 3
}

TEST(GeneEncoding, EmptyGenes) {
    baysor::MoleculeData data;
    std::vector<std::string> genes = {"A", "", "B"};
    baysor::encode_genes(data, genes);

    // Empty string is skipped from gene_names
    ASSERT_EQ(data.gene_names.size(), 2u);
    EXPECT_EQ(data.gene_names[0], "A");
    EXPECT_EQ(data.gene_names[1], "B");

    // Empty string maps to 0
    EXPECT_EQ(data.gene[0], 1);
    EXPECT_EQ(data.gene[1], 0);
    EXPECT_EQ(data.gene[2], 2);
}

// ============================================================================
// Gene filtering by count
// ============================================================================

TEST(GeneFiltering, FilterByCount) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3, 4, 5, 6, 7};
    data.y = {1, 2, 3, 4, 5, 6, 7};
    std::vector<std::string> genes = {"A", "A", "A", "B", "C", "C", "C"};
    baysor::encode_genes(data, genes);

    // A: 3, B: 1, C: 3 — filter with min=2 should remove B
    baysor::filter_genes_by_count(data, 2);

    EXPECT_EQ(data.n_molecules(), 6);
    EXPECT_EQ(data.n_genes(), 2);
    EXPECT_EQ(data.gene_names[0], "A");
    EXPECT_EQ(data.gene_names[1], "C");
}

TEST(GeneFiltering, FilterByCountNoOp) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3};
    data.y = {1, 2, 3};
    std::vector<std::string> genes = {"A", "B", "A"};
    baysor::encode_genes(data, genes);

    // A: 2, B: 1 — min=1 should not filter anything
    baysor::filter_genes_by_count(data, 1);
    EXPECT_EQ(data.n_molecules(), 3);
    EXPECT_EQ(data.n_genes(), 2);
}

TEST(GeneFiltering, FilterByCountCompactsMetadata) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3, 4};
    data.y = {1, 2, 3, 4};
    data.confidence = {0.9, 0.8, 0.7, 0.6};
    std::vector<std::string> genes = {"A", "B", "A", "A"};
    baysor::encode_genes(data, genes);

    // A: 3, B: 1 — filter min=2 removes molecule 1 (index 1, gene B)
    baysor::filter_genes_by_count(data, 2);

    EXPECT_EQ(data.n_molecules(), 3);
    ASSERT_EQ(data.confidence.size(), 3u);
    EXPECT_DOUBLE_EQ(data.confidence[0], 0.9);
    EXPECT_DOUBLE_EQ(data.confidence[1], 0.7);
    EXPECT_DOUBLE_EQ(data.confidence[2], 0.6);
}

// ============================================================================
// Gene filtering by pattern
// ============================================================================

TEST(GeneFiltering, FilterByPattern) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3, 4, 5};
    data.y = {1, 2, 3, 4, 5};
    std::vector<std::string> genes = {"Blank1", "GeneA", "Blank2", "GeneB", "MALAT1"};
    baysor::encode_genes(data, genes);

    // Filter genes matching "Blank*" pattern
    baysor::filter_genes_by_pattern(data, {"Blank*"});

    EXPECT_EQ(data.n_molecules(), 3);
    // Remaining genes: GeneA, GeneB, MALAT1
    EXPECT_EQ(data.n_genes(), 3);
}

TEST(GeneFiltering, FilterByMultiplePatterns) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3, 4, 5};
    data.y = {1, 2, 3, 4, 5};
    std::vector<std::string> genes = {"Blank1", "GeneA", "Blank2", "GeneB", "MALAT1"};
    baysor::encode_genes(data, genes);

    baysor::filter_genes_by_pattern(data, {"Blank*", "MALAT1"});

    EXPECT_EQ(data.n_molecules(), 2);
    EXPECT_EQ(data.n_genes(), 2);
}

// ============================================================================
// Position matrix
// ============================================================================

TEST(MoleculeData, PositionMatrix2D) {
    baysor::MoleculeData data;
    data.x = {1.0, 2.0, 3.0};
    data.y = {4.0, 5.0, 6.0};

    auto mat = data.position_matrix();
    EXPECT_EQ(mat.rows(), 2);
    EXPECT_EQ(mat.cols(), 3);
    EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(mat(1, 2), 6.0);
}

TEST(MoleculeData, PositionMatrix3D) {
    baysor::MoleculeData data;
    data.x = {1.0, 2.0};
    data.y = {3.0, 4.0};
    data.z = {5.0, 6.0};

    auto mat = data.position_matrix();
    EXPECT_EQ(mat.rows(), 3);
    EXPECT_EQ(mat.cols(), 2);
    EXPECT_DOUBLE_EQ(mat(2, 1), 6.0);
}

TEST(MoleculeData, Dimensionality) {
    baysor::MoleculeData data;
    data.x = {1.0};
    data.y = {2.0};
    EXPECT_FALSE(data.is_3d());
    EXPECT_EQ(data.n_dims(), 2);

    data.z = {3.0};
    EXPECT_TRUE(data.is_3d());
    EXPECT_EQ(data.n_dims(), 3);
}

// ============================================================================
// CSV loading
// ============================================================================

TEST(DataLoading, LoadCSVBasic) {
    auto path = write_temp_csv(
        "x,y,gene\n"
        "1.0,2.0,A\n"
        "3.0,4.0,B\n"
        "5.0,6.0,A\n"
    );

    baysor::DataOptions opts;
    auto data = baysor::load_molecules(path, opts);

    EXPECT_EQ(data.n_molecules(), 3);
    EXPECT_EQ(data.n_genes(), 2);
    EXPECT_DOUBLE_EQ(data.x[0], 1.0);
    EXPECT_DOUBLE_EQ(data.y[2], 6.0);
    EXPECT_FALSE(data.is_3d());

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVWith3D) {
    auto path = write_temp_csv(
        "x,y,z,gene\n"
        "1.0,2.0,10.0,A\n"
        "3.0,4.0,20.0,B\n"
    );

    baysor::DataOptions opts;
    auto data = baysor::load_molecules(path, opts);

    EXPECT_TRUE(data.is_3d());
    EXPECT_EQ(data.n_dims(), 3);
    EXPECT_DOUBLE_EQ(data.z[1], 20.0);

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVDropConstantZ) {
    auto path = write_temp_csv(
        "x,y,z,gene\n"
        "1.0,2.0,0.0,A\n"
        "3.0,4.0,0.0,B\n"
    );

    baysor::DataOptions opts;
    auto data = baysor::load_molecules(path, opts);

    // z is constant -> should be dropped
    EXPECT_FALSE(data.is_3d());
    EXPECT_TRUE(data.z.empty());

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVForce2D) {
    auto path = write_temp_csv(
        "x,y,z,gene\n"
        "1.0,2.0,10.0,A\n"
        "3.0,4.0,20.0,B\n"
    );

    baysor::DataOptions opts;
    opts.force_2d = true;
    auto data = baysor::load_molecules(path, opts);

    EXPECT_FALSE(data.is_3d());
    EXPECT_TRUE(data.z.empty());

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVWithConfidence) {
    auto path = write_temp_csv(
        "x,y,gene,confidence\n"
        "1.0,2.0,A,0.95\n"
        "3.0,4.0,B,0.80\n"
    );

    baysor::DataOptions opts;
    auto data = baysor::load_molecules(path, opts);

    ASSERT_EQ(data.confidence.size(), 2u);
    EXPECT_DOUBLE_EQ(data.confidence[0], 0.95);
    EXPECT_DOUBLE_EQ(data.confidence[1], 0.80);

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVCustomColumns) {
    auto path = write_temp_csv(
        "pos_x,pos_y,transcript\n"
        "1.0,2.0,GeneA\n"
        "3.0,4.0,GeneB\n"
    );

    baysor::DataOptions opts;
    opts.x_col = "pos_x";
    opts.y_col = "pos_y";
    opts.gene_col = "transcript";
    auto data = baysor::load_molecules(path, opts);

    EXPECT_EQ(data.n_molecules(), 2);
    EXPECT_DOUBLE_EQ(data.x[0], 1.0);
    EXPECT_EQ(data.gene_names[0], "GeneA");

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVWithGeneFilter) {
    auto path = write_temp_csv(
        "x,y,gene\n"
        "1.0,1.0,A\n"
        "2.0,2.0,A\n"
        "3.0,3.0,A\n"
        "4.0,4.0,B\n"
        "5.0,5.0,C\n"
        "6.0,6.0,C\n"
        "7.0,7.0,C\n"
    );

    baysor::DataOptions opts;
    opts.min_molecules_per_gene = 3;
    auto data = baysor::load_molecules(path, opts);

    // B has 1, C has 3, A has 3 — B filtered out
    EXPECT_EQ(data.n_molecules(), 6);
    EXPECT_EQ(data.n_genes(), 2);

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVWithExcludeGenes) {
    auto path = write_temp_csv(
        "x,y,gene\n"
        "1.0,1.0,Blank1\n"
        "2.0,2.0,GeneA\n"
        "3.0,3.0,Blank2\n"
        "4.0,4.0,GeneB\n"
    );

    baysor::DataOptions opts;
    opts.exclude_genes = "Blank*";
    auto data = baysor::load_molecules(path, opts);

    EXPECT_EQ(data.n_molecules(), 2);
    EXPECT_EQ(data.n_genes(), 2);

    std::remove(path.c_str());
}

// ============================================================================
// Prior segmentation
// ============================================================================

TEST(PriorSegmentation, DetectType) {
    EXPECT_EQ(baysor::detect_prior_seg_type(""), baysor::PriorSegType::None);
    EXPECT_EQ(baysor::detect_prior_seg_type(":cell_id"), baysor::PriorSegType::Column);
    EXPECT_EQ(baysor::detect_prior_seg_type("/path/to/mask.tiff"), baysor::PriorSegType::Image);
}

TEST(PriorSegmentation, FilterLabels) {
    std::vector<int> labels = {1, 1, 1, 2, 3, 3, 3, 3};
    baysor::filter_segmentation_labels(labels, 3);

    // Segment 1: 3 molecules (passes), segment 2: 1 (filtered), segment 3: 4 (passes)
    EXPECT_EQ(labels[0], 1);
    EXPECT_EQ(labels[3], 0);  // was 2, now filtered
    EXPECT_EQ(labels[4], 3);
}

TEST(PriorSegmentation, FilterLabelsNoOp) {
    std::vector<int> labels = {1, 1, 2, 2};
    baysor::filter_segmentation_labels(labels, 2);

    // Both segments have 2 molecules, should pass
    EXPECT_EQ(labels[0], 1);
    EXPECT_EQ(labels[2], 2);
}

TEST(PriorSegmentation, ParseFromColumn) {
    auto path = write_temp_csv(
        "x,y,gene,cell_id\n"
        "1.0,1.0,A,cell1\n"
        "2.0,2.0,A,cell1\n"
        "3.0,3.0,A,cell1\n"
        "4.0,4.0,B,cell2\n"
        "5.0,5.0,B,0\n"
        "6.0,6.0,B,cell2\n"
    );

    auto labels = baysor::parse_prior_from_column(path, "cell_id", "0", 0);

    ASSERT_EQ(labels.size(), 6u);
    // cell1 -> 1, cell2 -> 2 (alphabetically sorted), "0" -> 0 (unassigned)
    EXPECT_EQ(labels[0], 1);  // cell1
    EXPECT_EQ(labels[1], 1);  // cell1
    EXPECT_EQ(labels[2], 1);  // cell1
    EXPECT_EQ(labels[3], 2);  // cell2
    EXPECT_EQ(labels[4], 0);  // unassigned
    EXPECT_EQ(labels[5], 2);  // cell2

    std::remove(path.c_str());
}

TEST(PriorSegmentation, ParseFromColumnWithFilter) {
    auto path = write_temp_csv(
        "x,y,gene,cell_id\n"
        "1.0,1.0,A,cell1\n"
        "2.0,2.0,A,cell1\n"
        "3.0,3.0,A,cell1\n"
        "4.0,4.0,B,cell2\n"
        "5.0,5.0,B,0\n"
    );

    // cell2 has only 1 molecule, filter with min=2
    auto labels = baysor::parse_prior_from_column(path, "cell_id", "0", 2);

    EXPECT_EQ(labels[3], 0);  // cell2 filtered out

    std::remove(path.c_str());
}

TEST(PriorSegmentation, EstimateScale) {
    // Create a simple grid of 4 cell centers spaced 10 units apart
    // Cells at (0,0), (10,0), (0,10), (10,10)
    // Each cell has 5 molecules clustered around its center
    Eigen::MatrixXd pos(2, 20);
    std::vector<int> assignment(20);

    // Cell 1: around (0, 0)
    for (int i = 0; i < 5; ++i) {
        pos(0, i) = 0.0 + i * 0.1;
        pos(1, i) = 0.0 + i * 0.1;
        assignment[i] = 1;
    }
    // Cell 2: around (10, 0)
    for (int i = 0; i < 5; ++i) {
        pos(0, 5 + i) = 10.0 + i * 0.1;
        pos(1, 5 + i) = 0.0 + i * 0.1;
        assignment[5 + i] = 2;
    }
    // Cell 3: around (0, 10)
    for (int i = 0; i < 5; ++i) {
        pos(0, 10 + i) = 0.0 + i * 0.1;
        pos(1, 10 + i) = 10.0 + i * 0.1;
        assignment[10 + i] = 3;
    }
    // Cell 4: around (10, 10)
    for (int i = 0; i < 5; ++i) {
        pos(0, 15 + i) = 10.0 + i * 0.1;
        pos(1, 15 + i) = 10.0 + i * 0.1;
        assignment[15 + i] = 4;
    }

    auto [scale, scale_std] = baysor::estimate_scale_from_assignment(pos, assignment, 3);

    // NN distance between centers ~10.0, so radius ~5.0
    EXPECT_NEAR(scale, 5.0, 0.5);
    EXPECT_GE(scale_std, 0.0);
}

// ============================================================================
// Top-level prior segmentation loading
// ============================================================================

TEST(PriorSegmentation, LoadFromColumn) {
    auto path = write_temp_csv(
        "x,y,gene,cell_id\n"
        "1.0,1.0,A,cell1\n"
        "2.0,2.0,A,cell1\n"
        "3.0,3.0,A,cell1\n"
        "10.0,10.0,B,cell2\n"
        "11.0,11.0,B,cell2\n"
        "12.0,12.0,B,cell2\n"
        "20.0,20.0,A,cell3\n"
        "21.0,21.0,A,cell3\n"
        "22.0,22.0,A,cell3\n"
    );

    baysor::DataOptions opts;
    auto data = baysor::load_molecules(path, opts);

    auto [scale, scale_std] = baysor::load_prior_segmentation(
        data, ":cell_id", path, "0", 0, 3, true);

    ASSERT_EQ(data.prior_segmentation.size(), 9u);
    EXPECT_GT(data.prior_segmentation[0], 0);  // assigned
    EXPECT_GT(scale, 0.0);

    std::remove(path.c_str());
}

TEST(PriorSegmentation, LoadNone) {
    baysor::MoleculeData data;
    data.x = {1.0};
    data.y = {1.0};

    auto [scale, scale_std] = baysor::load_prior_segmentation(
        data, "", "dummy", "0", 0, 3, false);

    EXPECT_DOUBLE_EQ(scale, -1.0);
    EXPECT_TRUE(data.prior_segmentation.empty());
}

// ============================================================================
// Options
// ============================================================================

TEST(Options, ParseScaleStd) {
    EXPECT_DOUBLE_EQ(baysor::parse_scale_std("25%", 10.0), 2.5);
    EXPECT_DOUBLE_EQ(baysor::parse_scale_std("50%", 10.0), 5.0);
    EXPECT_DOUBLE_EQ(baysor::parse_scale_std("3.0", 10.0), 3.0);
    EXPECT_DOUBLE_EQ(baysor::parse_scale_std("", 10.0), 2.5);  // default 25%
}

TEST(Options, DefaultParamValues) {
    EXPECT_EQ(baysor::default_param_value("min_molecules_per_segment", 12), 3);
    EXPECT_EQ(baysor::default_param_value("confidence_nn_id", 10), 6);
    EXPECT_EQ(baysor::default_param_value("n_cells_init", 10, 1000), 200);
}

TEST(Options, FillDataOptions) {
    baysor::DataOptions opts;
    opts.min_molecules_per_cell = 20;
    baysor::fill_and_check_data_options(opts);

    EXPECT_GT(opts.min_molecules_per_segment, 0);
    EXPECT_GT(opts.confidence_nn_id, 0);
}

TEST(Options, LoadConfigFromToml) {
    auto path = write_temp_csv(
        "[data]\n"
        "x = \"pos_x\"\n"
        "min_molecules_per_cell = 30\n"
        "\n"
        "[segmentation]\n"
        "scale = 7.5\n"
        "n_clusters = 6\n",
        ".toml"
    );

    auto opts = baysor::load_config(path);

    EXPECT_EQ(opts.data.x_col, "pos_x");
    EXPECT_EQ(opts.data.min_molecules_per_cell, 30);
    EXPECT_DOUBLE_EQ(opts.segmentation.scale, 7.5);
    EXPECT_EQ(opts.segmentation.n_clusters, 6);

    std::remove(path.c_str());
}

// ============================================================================
// Utility functions
// ============================================================================

TEST(Utils, CountArray) {
    std::vector<int> values = {1, 2, 3, 1, 2, 1};
    auto counts = baysor::count_array(values);

    ASSERT_EQ(counts.size(), 3u);
    EXPECT_EQ(counts[0], 3);  // value 1
    EXPECT_EQ(counts[1], 2);  // value 2
    EXPECT_EQ(counts[2], 1);  // value 3
}

TEST(Utils, SplitIds) {
    std::vector<int> factor = {2, 1, 3, 1, 2};
    auto groups = baysor::split_ids(factor);

    ASSERT_EQ(groups.size(), 3u);
    // Group 1: indices where factor == 1
    ASSERT_EQ(groups[0].size(), 2u);
    EXPECT_EQ(groups[0][0], 1);
    EXPECT_EQ(groups[0][1], 3);
    // Group 2: indices where factor == 2
    ASSERT_EQ(groups[1].size(), 2u);
    EXPECT_EQ(groups[1][0], 0);
    EXPECT_EQ(groups[1][1], 4);
}

TEST(Utils, SplitStringList) {
    auto result = baysor::split_string_list("Blank*, MALAT1 , NegCtrl");
    ASSERT_EQ(result.size(), 3u);
    EXPECT_EQ(result[0], "Blank*");
    EXPECT_EQ(result[1], "MALAT1");
    EXPECT_EQ(result[2], "NegCtrl");
}

TEST(Utils, Strip) {
    EXPECT_EQ(baysor::strip("  hello  "), "hello");
    EXPECT_EQ(baysor::strip("\thello\n"), "hello");
    EXPECT_EQ(baysor::strip(""), "");
}

// ============================================================================
// KNN (nanoflann-based)
// ============================================================================

TEST(KNN, SmallGrid) {
    // 4 points on a unit square
    Eigen::MatrixXd pts(2, 4);
    pts.col(0) << 0.0, 0.0;
    pts.col(1) << 1.0, 0.0;
    pts.col(2) << 0.0, 1.0;
    pts.col(3) << 1.0, 1.0;

    auto result = baysor::knn_parallel(pts, pts, 3, true);

    ASSERT_EQ(result.indices.size(), 4u);
    ASSERT_EQ(result.distances.size(), 4u);

    // For point (0,0), nearest is itself (dist=0), then (1,0) or (0,1) (dist=1)
    EXPECT_EQ(result.indices[0][0], 0);  // self
    EXPECT_NEAR(result.distances[0][0], 0.0, 1e-10);
    EXPECT_NEAR(result.distances[0][1], 1.0, 1e-10);  // adjacent
    EXPECT_NEAR(result.distances[0][2], 1.0, 1e-10);  // adjacent
}

TEST(KNN, SelfQuery) {
    Eigen::MatrixXd pts(2, 3);
    pts.col(0) << 0.0, 0.0;
    pts.col(1) << 3.0, 0.0;
    pts.col(2) << 0.0, 4.0;

    auto result = baysor::knn_parallel(pts, pts, 2, true);

    // First neighbor is always self
    for (int i = 0; i < 3; ++i) {
        EXPECT_EQ(result.indices[i][0], i);
        EXPECT_NEAR(result.distances[i][0], 0.0, 1e-10);
    }
}

TEST(KNN, 3D) {
    Eigen::MatrixXd pts(3, 2);
    pts.col(0) << 0.0, 0.0, 0.0;
    pts.col(1) << 3.0, 4.0, 0.0;

    auto result = baysor::knn_parallel(pts, pts, 2, true);

    EXPECT_NEAR(result.distances[0][1], 5.0, 1e-10);
    EXPECT_NEAR(result.distances[1][1], 5.0, 1e-10);
}

// ============================================================================
// count_array_sparse
// ============================================================================

TEST(SparseCount, Basic) {
    // Gene IDs (1-based): {1, 2, 1, 3, 2}
    int values[] = {1, 2, 1, 3, 2};
    auto sv = baysor::count_array_sparse(values, 5, 3);

    EXPECT_EQ(sv.size(), 3);
    EXPECT_FLOAT_EQ(sv.coeff(0), 2.0f);  // gene 1: 2 counts
    EXPECT_FLOAT_EQ(sv.coeff(1), 2.0f);  // gene 2: 2 counts
    EXPECT_FLOAT_EQ(sv.coeff(2), 1.0f);  // gene 3: 1 count
}

TEST(SparseCount, Weighted) {
    int values[] = {1, 1, 2};
    double weights[] = {0.5, 1.5, 3.0};
    auto sv = baysor::count_array_sparse(values, 3, 2, weights);

    EXPECT_FLOAT_EQ(sv.coeff(0), 2.0f);  // gene 1: 0.5 + 1.5
    EXPECT_FLOAT_EQ(sv.coeff(1), 3.0f);  // gene 2: 3.0
}

TEST(SparseCount, Normalized) {
    int values[] = {1, 2, 2};
    auto sv = baysor::count_array_sparse(values, 3, 2, nullptr, true);

    // gene 1: 1/3, gene 2: 2/3
    EXPECT_NEAR(sv.coeff(0), 1.0f / 3.0f, 1e-5f);
    EXPECT_NEAR(sv.coeff(1), 2.0f / 3.0f, 1e-5f);
}

// ============================================================================
// Triangulation
// ============================================================================

TEST(Triangulation, NormalizePoints) {
    Eigen::MatrixXd pts(2, 4);
    pts.col(0) << 0.0, 0.0;
    pts.col(1) << 100.0, 0.0;
    pts.col(2) << 0.0, 100.0;
    pts.col(3) << 100.0, 100.0;

    auto norm = baysor::normalize_points(pts);

    // All values should be in [1.0, 2.0]
    EXPECT_GE(norm.minCoeff(), 1.0);
    EXPECT_LE(norm.maxCoeff(), 2.0);
}

TEST(Triangulation, AdjacencyList2D) {
    // Create a small 2D grid: 4 points on a unit square
    Eigen::MatrixXd pts(2, 4);
    pts.col(0) << 0.0, 0.0;
    pts.col(1) << 1.0, 0.0;
    pts.col(2) << 0.0, 1.0;
    pts.col(3) << 1.0, 1.0;

    auto result = baysor::adjacency_list(pts, /*filter=*/false);

    // Delaunay of 4 points on a square should produce 5 edges
    // (4 sides + 1 diagonal)
    EXPECT_GE(static_cast<int>(result.edge_src.size()), 4);
    EXPECT_LE(static_cast<int>(result.edge_src.size()), 6);

    // All edges should have positive distance
    for (double d : result.edge_dists) {
        EXPECT_GT(d, 0.0);
    }
}

TEST(Triangulation, AdjacencyList3DKNN) {
    // 3D should use KNN
    Eigen::MatrixXd pts(3, 5);
    pts.col(0) << 0.0, 0.0, 0.0;
    pts.col(1) << 1.0, 0.0, 0.0;
    pts.col(2) << 0.0, 1.0, 0.0;
    pts.col(3) << 0.0, 0.0, 1.0;
    pts.col(4) << 1.0, 1.0, 1.0;

    auto result = baysor::adjacency_list(pts, /*filter=*/false);

    EXPECT_GT(static_cast<int>(result.edge_src.size()), 0);
}

// ============================================================================
// AdjList::from_edge_list
// ============================================================================

TEST(AdjList, FromEdgeList) {
    // 3 vertices, 2 edges: 0-1, 1-2
    int src[] = {0, 1};
    int dst[] = {1, 2};
    double wts[] = {0.5, 0.8};

    auto adj = baysor::AdjList::from_edge_list(src, dst, wts, 2, 3);

    EXPECT_EQ(adj.n_molecules(), 3);

    // Vertex 0: neighbors should include 1
    EXPECT_GE(adj.neighbor_count(0), 1);
    // Vertex 1: neighbors should include 0 and 2
    EXPECT_GE(adj.neighbor_count(1), 2);
    // Vertex 2: neighbors should include 1
    EXPECT_GE(adj.neighbor_count(2), 1);
}

TEST(AdjList, EmptyGraph) {
    auto adj = baysor::AdjList::from_edge_list(nullptr, nullptr, nullptr, 0, 5);

    EXPECT_EQ(adj.n_molecules(), 5);
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(adj.neighbor_count(i), 0);
    }
}

// ============================================================================
// build_molecule_graph
// ============================================================================

TEST(MoleculeGraph, Basic) {
    // Use a 2D grid to avoid degenerate collinear triangulation
    baysor::MoleculeData data;
    // 3x3 grid
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            data.x.push_back(i * 1.0);
            data.y.push_back(j * 1.0);
            data.gene.push_back(1);
        }
    }

    auto adj = baysor::build_molecule_graph(data);

    EXPECT_EQ(adj.n_molecules(), 9);
    // Center point (index 4) should have neighbors
    EXPECT_GE(adj.neighbor_count(4), 2);
    // Total edges should be reasonable
    EXPECT_GT(adj.nnz(), 0);
}

// ============================================================================
// Noise estimation
// ============================================================================

TEST(NoiseEstimation, FitBimodal) {
    // Create synthetic bimodal data: signal (low distances) + noise (high distances)
    std::vector<double> edge_lengths;
    std::mt19937 rng(42);

    // Signal: Normal(2.0, 0.5) — 100 points
    std::normal_distribution<double> signal_dist(2.0, 0.5);
    for (int i = 0; i < 100; ++i) {
        edge_lengths.push_back(std::max(0.1, signal_dist(rng)));
    }
    // Noise: Normal(8.0, 1.0) — 50 points
    std::normal_distribution<double> noise_dist(8.0, 1.0);
    for (int i = 0; i < 50; ++i) {
        edge_lengths.push_back(noise_dist(rng));
    }

    int n = static_cast<int>(edge_lengths.size());

    // Build a simple fully-connected adjacency for testing
    // (each point connected to all others — won't happen in practice but tests the EM)
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n && j < i + 5; ++j) {
            src.push_back(i);
            dst.push_back(j);
            wts.push_back(1.0);
        }
    }
    auto adj = baysor::AdjList::from_edge_list(src.data(), dst.data(), wts.data(),
                                                static_cast<int>(src.size()), n);

    auto result = baysor::fit_noise_probabilities(edge_lengths, adj, nullptr, 100, 0.005, false);

    // Signal mu should be close to 2.0, noise mu close to 8.0
    EXPECT_NEAR(result.signal_mu, 2.0, 1.5);
    EXPECT_NEAR(result.noise_mu, 8.0, 2.0);

    // Assignment probs should be n x 2
    EXPECT_EQ(result.assignment_probs.rows(), n);
    EXPECT_EQ(result.assignment_probs.cols(), 2);

    // Signal points should have high signal probability
    double avg_signal_prob = 0.0;
    for (int i = 0; i < 100; ++i) {
        avg_signal_prob += result.assignment_probs(i, 0);
    }
    avg_signal_prob /= 100.0;
    EXPECT_GT(avg_signal_prob, 0.5);

    // Noise points should have low signal probability
    double avg_noise_prob = 0.0;
    for (int i = 100; i < n; ++i) {
        avg_noise_prob += result.assignment_probs(i, 0);
    }
    avg_noise_prob /= 50.0;
    EXPECT_LT(avg_noise_prob, 0.5);
}

TEST(NoiseEstimation, AppendConfidence) {
    // Create a grid of molecules with varied spacing
    baysor::MoleculeData data;
    // Dense cluster (signal-like)
    for (int i = 0; i < 20; ++i) {
        data.x.push_back(i * 0.5);
        data.y.push_back(0.0);
        data.gene.push_back(1);
    }
    // Sparse points (noise-like)
    for (int i = 0; i < 10; ++i) {
        data.x.push_back(100.0 + i * 10.0);
        data.y.push_back(100.0);
        data.gene.push_back(1);
    }
    data.gene_names = {"GeneA"};

    baysor::append_confidence(data, 3);

    ASSERT_EQ(data.confidence.size(), 30u);

    // All confidence values should be in [0, 1]
    for (double c : data.confidence) {
        EXPECT_GE(c, 0.0);
        EXPECT_LE(c, 1.0);
    }

    // Dense cluster should generally have higher confidence than sparse points
    double avg_dense = 0.0;
    for (int i = 0; i < 20; ++i) avg_dense += data.confidence[i];
    avg_dense /= 20.0;

    double avg_sparse = 0.0;
    for (int i = 20; i < 30; ++i) avg_sparse += data.confidence[i];
    avg_sparse /= 10.0;

    EXPECT_GT(avg_dense, avg_sparse);
}

// ============================================================================
// Neighborhood composition
// ============================================================================

TEST(NeighborhoodComposition, CountMatrix) {
    // 6 molecules, 3 genes, on a line
    Eigen::MatrixXd pos(2, 6);
    for (int i = 0; i < 6; ++i) {
        pos(0, i) = i * 1.0;
        pos(1, i) = 0.0;
    }
    std::vector<int> genes = {1, 2, 3, 1, 2, 3};

    auto cm = baysor::neighborhood_count_matrix(pos, genes, 3, 3);

    // Should be 3 (genes) x 6 (molecules)
    EXPECT_EQ(cm.rows(), 3);
    EXPECT_EQ(cm.cols(), 6);

    // Each column should have non-zero entries (neighbors exist)
    for (int i = 0; i < 6; ++i) {
        double col_sum = 0;
        for (int j = 0; j < 3; ++j) {
            col_sum += cm.coeff(j, i);
        }
        EXPECT_GT(col_sum, 0.0);
    }
}

TEST(NeighborhoodComposition, EstimateGeneVectors) {
    // Create a small count matrix (3 genes x 10 molecules)
    Eigen::SparseMatrix<float> cm(3, 10);
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    for (int j = 0; j < 10; ++j) {
        for (int i = 0; i < 3; ++i) {
            if (dist(rng) > 0.3f) {
                cm.insert(i, j) = dist(rng);
            }
        }
    }
    cm.makeCompressed();

    std::vector<int> gene_ids = {1, 2, 3, 1, 2, 3, 1, 2, 3, 1};

    auto gene_vecs = baysor::estimate_gene_vectors(cm, gene_ids, 2, "ri", false);

    // Should be 2 (components) x 3 (genes)
    EXPECT_EQ(gene_vecs.rows(), 2);
    EXPECT_EQ(gene_vecs.cols(), 3);
}

TEST(NeighborhoodComposition, EstimateGeneVectorsPerMolecule) {
    Eigen::SparseMatrix<float> cm(3, 10);
    for (int j = 0; j < 10; ++j) {
        cm.insert(j % 3, j) = 1.0f;
    }
    cm.makeCompressed();

    std::vector<int> gene_ids(10, 1);
    auto mol_vecs = baysor::estimate_gene_vectors(cm, gene_ids, 5, "ri", true);

    // Should be 5 (components) x 10 (molecules)
    EXPECT_EQ(mol_vecs.rows(), 5);
    EXPECT_EQ(mol_vecs.cols(), 10);
}

// ============================================================================
// Color embedding
// ============================================================================

TEST(ColorUtils, NormalizeEmbedding) {
    Eigen::MatrixXd emb(3, 100);
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 100; ++j) {
            emb(i, j) = dist(rng);
        }
    }

    baysor::normalize_embedding_to_lab_range(emb);

    // L channel should be in [~10, ~90]
    for (int j = 0; j < 100; ++j) {
        EXPECT_GE(emb(0, j), 0.0);
        EXPECT_LE(emb(0, j), 100.0);
    }
    // a, b channels should be in [-100, 100]
    for (int j = 0; j < 100; ++j) {
        EXPECT_GE(emb(1, j), -110.0);
        EXPECT_LE(emb(1, j), 110.0);
        EXPECT_GE(emb(2, j), -110.0);
        EXPECT_LE(emb(2, j), 110.0);
    }
}

TEST(ColorUtils, EmbeddingToHex) {
    // Middle gray in LAB: L=50, a=0, b=0
    Eigen::MatrixXd lab(3, 2);
    lab.col(0) << 50.0, 0.0, 0.0;  // neutral gray
    lab.col(1) << 90.0, 0.0, 0.0;  // light gray

    auto hexes = baysor::embedding_to_hex(lab);

    ASSERT_EQ(hexes.size(), 2u);
    // Should start with '#' and be 7 chars
    EXPECT_EQ(hexes[0].size(), 7u);
    EXPECT_EQ(hexes[0][0], '#');
    EXPECT_EQ(hexes[1].size(), 7u);
    EXPECT_EQ(hexes[1][0], '#');
}

TEST(ColorUtils, GeneCompositionColorEmbedding) {
    // Create synthetic mol_vecs: 10 components x 50 molecules
    Eigen::MatrixXd mol_vecs(10, 50);
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 50; ++j) {
            mol_vecs(i, j) = dist(rng);
        }
    }
    std::vector<double> confidence(50, 0.9);

    auto colors = baysor::gene_composition_color_embedding(mol_vecs, confidence, 30);

    ASSERT_EQ(colors.size(), 50u);
    for (const auto& c : colors) {
        EXPECT_EQ(c.size(), 7u);
        EXPECT_EQ(c[0], '#');
    }
}
