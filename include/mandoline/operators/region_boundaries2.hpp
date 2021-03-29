
#pragma once

#include "mandoline/mesh2.hpp"

namespace mandoline::operators {
std::set<int> region_boundaries(const CutCellMesh<2>& ccm);
std::map<int, bool> boundary_from_cells(const CutCellMesh<2>& ccm,
                                        const std::set<int>& cells);
std::map<int, bool> region_boundaries(const CutCellMesh<2>& ccm,
                                      std::set<int>& regions);
std::map<int, bool> region_boundaries(const CutCellMesh<2>& ccm, int region);

}  // namespace mandoline::operators
