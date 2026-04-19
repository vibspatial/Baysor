#pragma once

#include <string>

namespace baysor {

struct XeniumManifestContext {
    std::string manifest_path;
    std::string dataset_dir;
    std::string transcripts_path;
};

bool is_xenium_manifest_path(const std::string& path);
XeniumManifestContext load_xenium_manifest_context(const std::string& manifest_path);

} // namespace baysor
