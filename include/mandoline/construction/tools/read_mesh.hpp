#pragma once
#include <filesystem>
#include <balsa/eigen/types.hpp>

namespace mandoline::construction::tools {

std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i> read_mesh(
    const std::filesystem::path& path, bool preprocess = true);
}
