#pragma once
#include <filesystem>
#include <mtao/types.hpp>

namespace mandoline::construction::tools {

std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> read_mesh(
    const std::filesystem::path& path, bool preprocess = true);
}
