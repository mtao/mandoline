#pragma once
#include "mandoline/mesh3.hpp"


namespace mandoline::operators {
    mtao::VecXd cell_volumes(const CutCellMesh<3>& ccm) ;
    mtao::VecXd face_volumes(const CutCellMesh<3>& ccm, bool from_triangulation = false);
    mtao::VecXd dual_edge_lengths(const CutCellMesh<3>& ccm);
    mtao::VecXd dual_hodge2(const CutCellMesh<3>& ccm);
    mtao::VecXd primal_hodge2(const CutCellMesh<3>& ccm);
    mtao::VecXd dual_hodge3(const CutCellMesh<3>& ccm);
    mtao::VecXd primal_hodge3(const CutCellMesh<3>& ccm);
}
