#include "baysor/utils/xenium.h"

#include <fstream>
#include <stdexcept>

namespace baysor {

namespace {

std::string dirname_of(const std::string& path) {
    auto pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return ".";
    if (pos == 0) return path.substr(0, 1);
    return path.substr(0, pos);
}

std::string join_path(const std::string& a, const std::string& b) {
    if (b.empty()) return a;
    if (!b.empty() && (b[0] == '/' || b[0] == '\\')) return b;
    if (a.empty() || a == ".") return b;
    if (a.back() == '/' || a.back() == '\\') return a + b;
    return a + "/" + b;
}

bool file_exists(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    return static_cast<bool>(f);
}

} // namespace

bool is_xenium_manifest_path(const std::string& path) {
    constexpr const char* suffix = ".xenium";
    constexpr std::size_t suffix_len = 7;
    return path.size() >= suffix_len && path.substr(path.size() - suffix_len) == suffix;
}

XeniumManifestContext load_xenium_manifest_context(const std::string& manifest_path) {
    std::ifstream in(manifest_path);
    if (!in) {
        throw std::runtime_error("Could not open Xenium manifest '" + manifest_path + "'");
    }

    const std::string dataset_dir = dirname_of(manifest_path);

    std::string transcripts_path = join_path(dataset_dir, "transcripts.parquet");
    if (!file_exists(transcripts_path)) {
        transcripts_path = join_path(dataset_dir, "transcripts.csv.gz");
    }
    if (!file_exists(transcripts_path)) {
        throw std::runtime_error(
            "Could not locate transcripts.parquet or transcripts.csv.gz next to Xenium manifest '" +
            manifest_path + "'");
    }

    return XeniumManifestContext{
        manifest_path,
        dataset_dir,
        transcripts_path
    };
}

} // namespace baysor
