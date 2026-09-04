#include <gtest/gtest.h>

#include "baysor/segmentation/segmentation.h"
#include "baysor/utils/options.h"

#include <omp.h>

#include <atomic>
#include <cstdint>
#include <filesystem>
#include <future>
#include <string>
#include <thread>
#include <variant>
#include <vector>

#ifndef BAYSOR_TEST_FIXTURE_DIR
#error "BAYSOR_TEST_FIXTURE_DIR must identify the native baseline fixture"
#endif

namespace {

namespace fs = std::filesystem;

class ScopedTestOpenmpState {
public:
    ScopedTestOpenmpState()
        : max_threads_(omp_get_max_threads()), dynamic_(omp_get_dynamic()) {}

    ~ScopedTestOpenmpState() {
        omp_set_num_threads(max_threads_);
        omp_set_dynamic(dynamic_);
    }

    ScopedTestOpenmpState(const ScopedTestOpenmpState&) = delete;
    ScopedTestOpenmpState& operator=(const ScopedTestOpenmpState&) = delete;

private:
    int max_threads_;
    int dynamic_;
};

baysor::SegmentationRequest make_lifecycle_request() {
    const fs::path fixture_dir(BAYSOR_TEST_FIXTURE_DIR);
    const baysor::RunOptions config = baysor::load_config(
        (fixture_dir / "config.toml").string());

    baysor::SegmentationRequest request;
    request.molecules.path = (fixture_dir / "molecules.csv").string();
    request.molecules.options = config.molecules;
    request.prior = config.prior;
    request.segmentation = config.segmentation;
    request.neighborhood_composition.neighborhood_size =
        config.plotting.gene_composition_neighborhood;
    request.neighborhood_composition.method = config.plotting.ncv_method;
    request.requested_products = baysor::SegmentationProducts::none();
    request.requested_products.molecule_assignments = true;
    request.random_seed = baysor::kDefaultSegmentationSeed;
    request.execution.native_threads = 1;
    request.execution.use_arrow_threads = false;
    return request;
}

const baysor::SegmentationResult& require_result(
    const baysor::SegmentationOutcome& outcome
) {
    EXPECT_TRUE(std::holds_alternative<baysor::SegmentationResult>(outcome));
    return std::get<baysor::SegmentationResult>(outcome);
}

TEST(SegmentationLifecycle, InvalidAndInputFailuresHaveStructuredCategories) {
    baysor::SegmentationRequest invalid;
    try {
        (void)baysor::run_segmentation(invalid, baysor::CancellationToken{});
        FAIL() << "Expected an invalid-request error";
    } catch (const baysor::SegmentationError& error) {
        EXPECT_EQ(error.code(), baysor::SegmentationErrorCode::InvalidRequest);
    }

    auto missing_input = make_lifecycle_request();
    missing_input.molecules.path =
        (fs::path(BAYSOR_TEST_FIXTURE_DIR) / "missing-molecules.csv").string();
    try {
        (void)baysor::run_segmentation(missing_input, baysor::CancellationToken{});
        FAIL() << "Expected a molecule-input error";
    } catch (const baysor::SegmentationError& error) {
        EXPECT_EQ(error.code(), baysor::SegmentationErrorCode::MoleculeInput);
    }

    auto missing_prior = make_lifecycle_request();
    missing_prior.prior.type = baysor::PriorInputType::Boundary;
    missing_prior.prior.path =
        (fs::path(BAYSOR_TEST_FIXTURE_DIR) / "missing-boundaries.csv").string();
    try {
        (void)baysor::run_segmentation(missing_prior, baysor::CancellationToken{});
        FAIL() << "Expected a prior-input error";
    } catch (const baysor::SegmentationError& error) {
        EXPECT_EQ(error.code(), baysor::SegmentationErrorCode::PriorInput);
    }
}

TEST(SegmentationLifecycle, ZeroThreadsInheritsOpenmpConfigurationAndRecordsProvenance) {
    ScopedTestOpenmpState restore_openmp;
    omp_set_num_threads(2);
    const int inherited_max_threads = omp_get_max_threads();
    const bool inherited_dynamic = omp_get_dynamic() != 0;

    auto request = make_lifecycle_request();
    request.execution.native_threads = 0;
    const auto outcome = baysor::run_segmentation(request, baysor::CancellationToken{});
    const auto& result = require_result(outcome);

    EXPECT_EQ(result.resolved_options.execution.native_threads, 0);
    EXPECT_FALSE(result.resolved_options.execution.use_arrow_threads);
    EXPECT_EQ(result.provenance.requested_native_threads, 0);
    EXPECT_EQ(result.provenance.configured_openmp_max_threads, inherited_max_threads);
    EXPECT_EQ(result.provenance.openmp_dynamic_enabled, inherited_dynamic);
    EXPECT_FALSE(result.provenance.arrow_threads_enabled);
    EXPECT_EQ(omp_get_max_threads(), inherited_max_threads);
    EXPECT_EQ(omp_get_dynamic() != 0, inherited_dynamic);
}

TEST(SegmentationLifecycle, PositiveThreadRequestIsRestoredAfterSuccess) {
    ScopedTestOpenmpState restore_openmp;
    omp_set_num_threads(3);
    const int caller_max_threads = omp_get_max_threads();
    const bool caller_dynamic = omp_get_dynamic() != 0;

    auto request = make_lifecycle_request();
    request.execution.native_threads = 1;
    const auto outcome = baysor::run_segmentation(request, baysor::CancellationToken{});
    const auto& result = require_result(outcome);

    EXPECT_EQ(result.provenance.requested_native_threads, 1);
    EXPECT_EQ(result.provenance.configured_openmp_max_threads, 1);
    EXPECT_EQ(result.provenance.openmp_dynamic_enabled, caller_dynamic);
    EXPECT_EQ(omp_get_max_threads(), caller_max_threads);
    EXPECT_EQ(omp_get_dynamic() != 0, caller_dynamic);
}

TEST(SegmentationLifecycle, PositiveThreadRequestIsRestoredAfterException) {
    ScopedTestOpenmpState restore_openmp;
    omp_set_num_threads(3);
    const int caller_max_threads = omp_get_max_threads();
    const bool caller_dynamic = omp_get_dynamic() != 0;

    auto request = make_lifecycle_request();
    request.molecules.path =
        (fs::path(BAYSOR_TEST_FIXTURE_DIR) / "missing-molecules.csv").string();
    request.execution.native_threads = 1;
    EXPECT_THROW(
        (void)baysor::run_segmentation(request, baysor::CancellationToken{}),
        baysor::SegmentationError);

    EXPECT_EQ(omp_get_max_threads(), caller_max_threads);
    EXPECT_EQ(omp_get_dynamic() != 0, caller_dynamic);
}

TEST(SegmentationLifecycle, PositiveThreadRequestIsRestoredAfterCancellation) {
    ScopedTestOpenmpState restore_openmp;
    omp_set_num_threads(3);
    const int caller_max_threads = omp_get_max_threads();
    const bool caller_dynamic = omp_get_dynamic() != 0;

    auto request = make_lifecycle_request();
    request.execution.native_threads = 1;
    baysor::CancellationSource source;
    ASSERT_TRUE(source.request_cancellation());

    const auto outcome = baysor::run_segmentation(request, source.token());
    EXPECT_TRUE(std::holds_alternative<baysor::SegmentationCancelled>(outcome));
    EXPECT_FALSE(std::holds_alternative<baysor::SegmentationResult>(outcome));
    EXPECT_EQ(omp_get_max_threads(), caller_max_threads);
    EXPECT_EQ(omp_get_dynamic() != 0, caller_dynamic);
}

TEST(SegmentationLifecycle, ActiveOpenmpRegionIsRejectedBeforeConfigurationChanges) {
    ScopedTestOpenmpState restore_openmp;
    omp_set_dynamic(0);
    auto request = make_lifecycle_request();
    request.execution.native_threads = 1;
    std::atomic<bool> rejected{false};
    std::atomic<int> error_code{-1};

#pragma omp parallel num_threads(2)
    {
#pragma omp single
        {
            try {
                (void)baysor::run_segmentation(request, baysor::CancellationToken{});
            } catch (const baysor::SegmentationError& error) {
                rejected.store(true, std::memory_order_release);
                error_code.store(static_cast<int>(error.code()), std::memory_order_release);
            }
        }
    }

    EXPECT_TRUE(rejected.load(std::memory_order_acquire));
    EXPECT_EQ(
        error_code.load(std::memory_order_acquire),
        static_cast<int>(baysor::SegmentationErrorCode::UnsupportedExecutionContext));
}

TEST(SegmentationLifecycle, ConcurrentCallsFailFastInsteadOfOversubscribing) {
    auto request = make_lifecycle_request();
    request.segmentation.iters = 1000000;

    baysor::CancellationSource cancellation_source;
    std::promise<void> start_promise;
    const auto start = start_promise.get_future().share();
    std::atomic<int> cancelled{0};
    std::atomic<int> rejected{0};
    std::atomic<int> unexpected{0};

    const auto invoke = [&]() {
        start.wait();
        try {
            const auto outcome = baysor::run_segmentation(
                request, cancellation_source.token());
            if (std::holds_alternative<baysor::SegmentationCancelled>(outcome)) {
                cancelled.fetch_add(1, std::memory_order_relaxed);
            } else {
                unexpected.fetch_add(1, std::memory_order_relaxed);
            }
        } catch (const baysor::SegmentationError& error) {
            if (error.code() == baysor::SegmentationErrorCode::UnsupportedExecutionContext) {
                rejected.fetch_add(1, std::memory_order_relaxed);
                cancellation_source.request_cancellation();
            } else {
                unexpected.fetch_add(1, std::memory_order_relaxed);
                cancellation_source.request_cancellation();
            }
        } catch (...) {
            unexpected.fetch_add(1, std::memory_order_relaxed);
            cancellation_source.request_cancellation();
        }
    };

    std::thread first(invoke);
    std::thread second(invoke);
    start_promise.set_value();
    first.join();
    second.join();

    EXPECT_EQ(rejected.load(std::memory_order_relaxed), 1);
    EXPECT_EQ(cancelled.load(std::memory_order_relaxed), 1);
    EXPECT_EQ(unexpected.load(std::memory_order_relaxed), 0);
}

TEST(SegmentationLifecycle, FailedAndCancelledCallsDoNotAdvanceLaterRunState) {
    const auto run_assignments = [](std::uint64_t seed) {
        auto request = make_lifecycle_request();
        request.random_seed = seed;
        const auto outcome = baysor::run_segmentation(
            request, baysor::CancellationToken{});
        return require_result(outcome).cell_assignments;
    };

    const auto first = run_assignments(baysor::kDefaultSegmentationSeed);

    baysor::CancellationSource cancellation_source;
    ASSERT_TRUE(cancellation_source.request_cancellation());
    const auto cancelled = baysor::run_segmentation(
        make_lifecycle_request(), cancellation_source.token());
    ASSERT_TRUE(std::holds_alternative<baysor::SegmentationCancelled>(cancelled));

    auto missing_input = make_lifecycle_request();
    missing_input.molecules.path =
        (fs::path(BAYSOR_TEST_FIXTURE_DIR) / "missing-molecules.csv").string();
    EXPECT_THROW(
        (void)baysor::run_segmentation(missing_input, baysor::CancellationToken{}),
        baysor::SegmentationError);

    const auto second = run_assignments(baysor::kDefaultSegmentationSeed);
    EXPECT_EQ(second, first);
}

} // namespace
