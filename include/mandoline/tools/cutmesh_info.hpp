#pragma once
#include "mandoline/mesh3.hpp"

namespace mandoline::tools {
void print_file_info(const std::string& filename);
void print_all_info(const CutCellMesh<3>& ccm);
void print_general_info(const CutCellMesh<3>& ccm);
void print_face_info(const CutCellMesh<3>& ccm);
void print_region_info(const CutCellMesh<3>& ccm);
void print_cell_info(const CutCellMesh<3>& ccm);
} // namespace mandoline::tools
