#pragma once
#include <nlohmann/json.hpp>

#include "mandoline/mesh3.hpp"

namespace mandoline::tools {
void print_file_info(const std::string& filename,
                     nlohmann::json* json = nullptr);
void print_all_info(const CutCellMesh<3>& ccm, nlohmann::json* json = nullptr);
void print_general_info(const CutCellMesh<3>& ccm,
                        nlohmann::json* json = nullptr);
void print_face_info(const CutCellMesh<3>& ccm, nlohmann::json* json = nullptr);
void print_region_info(const CutCellMesh<3>& ccm,
                       nlohmann::json* json = nullptr);
void print_volume_info(const CutCellMesh<3>& ccm,
                       nlohmann::json* json = nullptr);
void print_cell_info(const CutCellMesh<3>& ccm, nlohmann::json* json = nullptr);
void print_boundary_info(const CutCellMesh<3>& ccm, nlohmann::json* json = nullptr);

// something somewhat specialized that was requested
void print_face_adjacent_grid_info(const CutCellMesh<3>& ccm,
                                   nlohmann::json* json = nullptr);
}  // namespace mandoline::tools
