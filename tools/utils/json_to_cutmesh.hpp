#pragma once 
#include <mandoline/mesh3.hpp>
#include <nlohmann/json.hpp>
#include <filesystem>


// if a filename is provided then a relative path is used
mandoline::CutCellMesh<3> json_to_cutmesh(const std::string& filename);
mandoline::CutCellMesh<3> json_to_cutmesh(const nlohmann::json& js, const std::filesystem::path& relative_path = {});
