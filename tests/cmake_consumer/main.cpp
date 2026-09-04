#include "baysor/segmentation/segmentation.h"

int main() {
    return baysor::derive_random_substream_seed(
        baysor::kDefaultSegmentationSeed,
        baysor::RandomSubstream::CoreScientific
    ) == baysor::kDefaultSegmentationSeed ? 0 : 1;
}
