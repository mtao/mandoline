
#pragma once

#include "mandoline/mesh3.hpp"

namespace mandoline::operators {
std::set<int> region_boundaries(const CutCellMesh<3>& ccm);
std::map<int, bool> boundary_from_cells(const CutCellMesh<3>& ccm,
                                        const std::set<int>& cells);
std::map<int, bool> region_boundaries(const CutCellMesh<3>& ccm,
                                      std::set<int>& regions);
std::map<int, bool> region_boundaries(const CutCellMesh<3>& ccm, int region);

}  // namespace mandoline::operators
