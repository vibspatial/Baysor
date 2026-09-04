#include <gtest/gtest.h>

#include "baysor/segmentation/segmentation.h"

#include <atomic>
#include <chrono>
#include <future>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>

namespace {

using SegmentationOperation = baysor::SegmentationOutcome (*)(
    const baysor::SegmentationRequest&,
    const baysor::CancellationToken&
);

static_assert(std::is_same_v<decltype(&baysor::run_segmentation), SegmentationOperation>);
static_assert(std::is_copy_constructible_v<baysor::CancellationToken>);
static_assert(std::is_move_constructible_v<baysor::SegmentationResult>);

TEST(SegmentationContract, RequestDefaultsDescribeACompletePathOrientedRun) {
    const baysor::SegmentationRequest request;

    EXPECT_TRUE(request.molecules.path.empty());
    EXPECT_EQ(request.molecules.options.x_col, "x");
    EXPECT_EQ(request.molecules.options.y_col, "y");
    EXPECT_EQ(request.molecules.options.gene_col, "gene");
    EXPECT_EQ(request.prior.type, baysor::PriorInputType::None);
    EXPECT_EQ(request.random_seed, baysor::kDefaultSegmentationSeed);
    EXPECT_EQ(request.execution.native_threads, 0);
    EXPECT_TRUE(request.execution.use_arrow_threads);

    EXPECT_TRUE(request.requested_products.molecule_assignments);
    EXPECT_TRUE(request.requested_products.molecule_confidence);
    EXPECT_TRUE(request.requested_products.assignment_confidence);
    EXPECT_TRUE(request.requested_products.molecule_clusters);
    EXPECT_TRUE(request.requested_products.neighborhood_composition_colors);
    EXPECT_TRUE(request.requested_products.cell_statistics);
    EXPECT_TRUE(request.requested_products.boundaries);
    EXPECT_TRUE(request.requested_products.count_matrix);
    EXPECT_TRUE(request.requested_products.diagnostics);
}

TEST(SegmentationContract, ResultDefaultsDoNotClaimUnproducedProducts) {
    const baysor::SegmentationResult result;

    EXPECT_FALSE(result.produced_products.molecule_assignments);
    EXPECT_FALSE(result.produced_products.molecule_confidence);
    EXPECT_FALSE(result.produced_products.assignment_confidence);
    EXPECT_FALSE(result.produced_products.molecule_clusters);
    EXPECT_FALSE(result.produced_products.neighborhood_composition_colors);
    EXPECT_FALSE(result.produced_products.cell_statistics);
    EXPECT_FALSE(result.produced_products.boundaries);
    EXPECT_FALSE(result.produced_products.count_matrix);
    EXPECT_FALSE(result.produced_products.diagnostics);
}

TEST(SegmentationContract, RandomSubstreamContractPreservesCompatibilityDefaults) {
    EXPECT_EQ(baysor::kRandomSubstreamContractVersion, 1U);
    EXPECT_EQ(
        baysor::derive_random_substream_seed(1, baysor::RandomSubstream::CoreScientific),
        1U
    );
    EXPECT_EQ(
        baysor::derive_random_substream_seed(1, baysor::RandomSubstream::MoleculeClustering),
        42U
    );
    EXPECT_EQ(
        baysor::derive_random_substream_seed(1, baysor::RandomSubstream::NeighborhoodComposition),
        42U
    );
    EXPECT_EQ(
        baysor::derive_random_substream_seed(1, baysor::RandomSubstream::Diagnostics),
        42U
    );

    EXPECT_EQ(
        baysor::derive_random_substream_seed(7, baysor::RandomSubstream::CoreScientific),
        7191089600892374487ULL
    );
    EXPECT_EQ(
        baysor::derive_random_substream_seed(7, baysor::RandomSubstream::MoleculeClustering),
        309689372594955804ULL
    );
    EXPECT_EQ(
        baysor::derive_random_substream_seed(7, baysor::RandomSubstream::NeighborhoodComposition),
        16616101746815609346ULL
    );
    EXPECT_EQ(
        baysor::derive_random_substream_seed(7, baysor::RandomSubstream::Diagnostics),
        10753165928301472203ULL
    );
    EXPECT_STREQ(
        baysor::random_substream_name(baysor::RandomSubstream::NeighborhoodComposition),
        "neighborhood_composition"
    );
}

TEST(SegmentationContract, CancellationIsVisibleAcrossThreads) {
    baysor::CancellationSource source;
    const baysor::CancellationToken token = source.token();
    std::promise<void> observer_started;
    auto observer_ready = observer_started.get_future();
    std::atomic<bool> observer_saw_cancellation{false};

    std::thread observer([token, &observer_started, &observer_saw_cancellation]() {
        observer_started.set_value();
        const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
        while (!token.is_cancellation_requested() &&
               std::chrono::steady_clock::now() < deadline) {
            std::this_thread::yield();
        }
        observer_saw_cancellation.store(
            token.is_cancellation_requested(),
            std::memory_order_release
        );
    });

    observer_ready.wait();
    EXPECT_TRUE(source.request_cancellation());
    EXPECT_FALSE(source.request_cancellation());
    observer.join();

    EXPECT_TRUE(source.is_cancellation_requested());
    EXPECT_TRUE(token.is_cancellation_requested());
    EXPECT_TRUE(observer_saw_cancellation.load(std::memory_order_acquire));
}

TEST(SegmentationContract, TokenKeepsCancellationStateAlive) {
    baysor::CancellationToken token;
    {
        baysor::CancellationSource source;
        token = source.token();
        EXPECT_TRUE(source.request_cancellation());
    }

    EXPECT_TRUE(token.is_cancellation_requested());
    EXPECT_FALSE(baysor::CancellationToken{}.is_cancellation_requested());
}

TEST(SegmentationContract, ResultRetainsOwnedValuesAfterMove) {
    baysor::SegmentationResult source;
    source.molecules.x = {1.25};
    source.molecules.y = {2.5};
    source.molecules.gene = {1};
    source.molecules.gene_names = {"GeneA"};
    source.molecules.source_transcript_id = {1234};
    source.cell_assignments = {1};
    source.cell_ids = {"cell_1"};

    baysor::CellStatistics statistics;
    statistics.cell_ids = source.cell_ids;
    statistics.columns = {"x", "y"};
    statistics.values.resize(1, 2);
    statistics.values << 1.25, 2.5;
    source.cell_statistics = std::move(statistics);

    baysor::PolygonCollection boundaries;
    Eigen::MatrixXd polygon(3, 2);
    polygon << 0.0, 0.0,
               1.0, 0.0,
               0.0, 1.0;
    boundaries.emplace("cell_1", std::move(polygon));
    source.boundaries_2d = std::move(boundaries);

    baysor::CellByGeneCounts counts;
    counts.cell_ids = source.cell_ids;
    counts.gene_names = source.molecules.gene_names;
    counts.values.resize(1, 1);
    counts.values.insert(0, 0) = 1.0F;
    source.count_matrix = std::move(counts);
    source.produced_products.molecule_assignments = true;
    source.produced_products.cell_statistics = true;
    source.produced_products.boundaries = true;
    source.produced_products.count_matrix = true;

    baysor::SegmentationResult result = std::move(source);

    ASSERT_EQ(result.molecules.source_transcript_id.size(), 1U);
    EXPECT_EQ(result.molecules.source_transcript_id.front(), 1234U);
    ASSERT_TRUE(result.cell_statistics.has_value());
    EXPECT_DOUBLE_EQ(result.cell_statistics->values(0, 0), 1.25);
    ASSERT_TRUE(result.boundaries_2d.has_value());
    EXPECT_DOUBLE_EQ(result.boundaries_2d->at("cell_1")(2, 1), 1.0);
    ASSERT_TRUE(result.count_matrix.has_value());
    EXPECT_FLOAT_EQ(result.count_matrix->values.coeff(0, 0), 1.0F);
}

TEST(SegmentationContract, CancelledOutcomeCannotBeMistakenForSuccess) {
    baysor::SegmentationOutcome cancelled = baysor::SegmentationCancelled{};
    EXPECT_TRUE(std::holds_alternative<baysor::SegmentationCancelled>(cancelled));
    EXPECT_FALSE(std::holds_alternative<baysor::SegmentationResult>(cancelled));

    baysor::SegmentationOutcome completed = baysor::SegmentationResult{};
    EXPECT_TRUE(std::holds_alternative<baysor::SegmentationResult>(completed));
    EXPECT_FALSE(std::holds_alternative<baysor::SegmentationCancelled>(completed));
}

TEST(SegmentationContract, NativeErrorsRetainTheirCategory) {
    const baysor::SegmentationError error(
        baysor::SegmentationErrorCode::InvalidRequest,
        "scale must be positive"
    );

    EXPECT_EQ(error.code(), baysor::SegmentationErrorCode::InvalidRequest);
    EXPECT_STREQ(error.what(), "scale must be positive");
}

} // namespace
