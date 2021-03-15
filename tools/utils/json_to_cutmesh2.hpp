#pragma once 
#include <mandoline/mesh2.hpp>
#include <nlohmann/json.hpp>
#include <filesystem>


// if a filename is provided then a relative path is used
mandoline::CutCellMesh<2> json_to_cutmesh2(const std::string& filename);
mandoline::CutCellMesh<2> json_to_cutmesh2(const nlohmann::json& js, const std::filesystem::path& relative_path = {});
