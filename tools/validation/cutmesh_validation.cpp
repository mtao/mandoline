#include "cutmesh_validation.hpp"

bool run_all(const mandoline::CutCellMesh<3>& ccm) {
    return pcwn_check(ccm) && faces_fully_utilized(ccm) &&
           grid_cells_fully_utilized(ccm) && equal_region_check(ccm) &&
           volume_check(ccm) && paired_boundary(ccm) &&
           exterior_cell_valence_counts(ccm);
}
