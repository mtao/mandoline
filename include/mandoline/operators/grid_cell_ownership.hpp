#pragma once
#include <set>
#include <vector>

#include "mandoline/mesh3.hpp"

namespace mandoline::operators {

std::vector<std::set<int>> grid_cell_ownership(const CutCellMesh<3>& ccm);
}
