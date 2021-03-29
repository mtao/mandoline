#pragma once 
#include "mandoline/mesh3.hpp"
#include <nlohmann/json.hpp>
#include <filesystem>

namespace mandoline::construction::tools {

// if a filename is provided then a relative path is used
CutCellMesh<3> json_to_cutmesh(const std::string& filename);
CutCellMesh<3> json_to_cutmesh(const nlohmann::json& js, const std::filesystem::path& relative_path = {});
}
