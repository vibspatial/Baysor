#include "baysor/segmentation/cancellation.h"

#include <atomic>
#include <utility>

namespace baysor {

namespace detail {

struct CancellationState {
    std::atomic<bool> requested{false};
};

} // namespace detail

CancellationToken::CancellationToken(std::shared_ptr<detail::CancellationState> state) noexcept
    : state_(std::move(state)) {}

bool CancellationToken::is_cancellation_requested() const noexcept {
    return state_ && state_->requested.load(std::memory_order_acquire);
}

CancellationSource::CancellationSource()
    : state_(std::make_shared<detail::CancellationState>()) {}

CancellationToken CancellationSource::token() const noexcept {
    return CancellationToken(state_);
}

bool CancellationSource::is_cancellation_requested() const noexcept {
    return state_ && state_->requested.load(std::memory_order_acquire);
}

bool CancellationSource::request_cancellation() noexcept {
    if (!state_) return false;

    bool expected = false;
    return state_->requested.compare_exchange_strong(
        expected,
        true,
        std::memory_order_acq_rel,
        std::memory_order_acquire
    );
}

} // namespace baysor
