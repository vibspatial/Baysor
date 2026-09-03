#include "baysor/segmentation/segmentation.h"

#include <utility>

namespace baysor {

namespace {

std::uint64_t splitmix64_finalizer(std::uint64_t value) noexcept {
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
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

} // namespace baysor
