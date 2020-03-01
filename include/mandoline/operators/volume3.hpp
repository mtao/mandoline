#pragma once
#include "mandoline/mesh3.hpp"


namespace mandoline::operators {
mtao::VecXd cell_volumes(const CutCellMesh<3> &ccm);
mtao::VecXd face_volumes(const CutCellMesh<3> &ccm, bool from_triangulation = false);
mtao::VecXd dual_edge_lengths(const CutCellMesh<3> &ccm);


// o^k/|a^k| = o^{n-k} / |a^{n-k}|
// H: o^k -> o^{n-k} = |a^{n-k}|/|a^k|


// dual_edge_lengths / face_volumes
// dual 1 -> primal 2
mtao::VecXd primal_hodge2(const CutCellMesh<3> &ccm);

// face_volumes / dual_edge_lengths
// primal 2 -> dual 1
mtao::VecXd dual_hodge2(const CutCellMesh<3> &ccm);

// 1.0 / cell_volumes
// primal 3 -> dual 0
mtao::VecXd primal_hodge3(const CutCellMesh<3> &ccm);

// cell_volumes
// dual 0 -> primal 3
mtao::VecXd dual_hodge3(const CutCellMesh<3> &ccm);
}// namespace mandoline::operators
