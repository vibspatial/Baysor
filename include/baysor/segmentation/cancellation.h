#pragma once

#include <memory>

namespace baysor {

namespace detail {
struct CancellationState;
}

class CancellationSource;

/// Read-only view of cooperative cancellation state.
///
/// Tokens are cheap to copy and may safely be observed from a different thread
/// than the source that requests cancellation. A default-constructed token has
/// no state and is never cancelled. Keeping a token alive also keeps its shared
/// state alive, even if every corresponding source has been destroyed.
class CancellationToken {
public:
    CancellationToken() noexcept = default;

    [[nodiscard]] bool is_cancellation_requested() const noexcept;

private:
    explicit CancellationToken(std::shared_ptr<detail::CancellationState> state) noexcept;

    std::shared_ptr<detail::CancellationState> state_;

    friend class CancellationSource;
};

/// Owner of the capability to request cooperative cancellation.
///
/// Copies share one cancellation state. request_cancellation() is idempotent and
/// returns true only for the call that first changes the state to cancelled.
class CancellationSource {
public:
    CancellationSource();

    [[nodiscard]] CancellationToken token() const noexcept;
    [[nodiscard]] bool is_cancellation_requested() const noexcept;
    bool request_cancellation() noexcept;

private:
    std::shared_ptr<detail::CancellationState> state_;
};

} // namespace baysor
