#pragma once
#include "mandoline/mesh3.hpp"


namespace mandoline::operators {
balsa::eigen::VecXd cell_volumes(const CutCellMesh<3> &ccm);
balsa::eigen::VecXd face_volumes(const CutCellMesh<3> &ccm, bool from_triangulation = false);
balsa::eigen::VecXd dual_edge_lengths(const CutCellMesh<3> &ccm);


// o^k/|a^k| = o^{n-k} / |a^{n-k}|
// H: o^k -> o^{n-k} = |a^{n-k}|/|a^k|


// dual_edge_lengths / face_volumes
// dual 1 -> primal 2
balsa::eigen::VecXd primal_hodge2(const CutCellMesh<3> &ccm);

// face_volumes / dual_edge_lengths
// primal 2 -> dual 1
balsa::eigen::VecXd dual_hodge2(const CutCellMesh<3> &ccm);

// 1.0 / cell_volumes
// primal 3 -> dual 0
balsa::eigen::VecXd primal_hodge3(const CutCellMesh<3> &ccm);

// cell_volumes
// dual 0 -> primal 3
balsa::eigen::VecXd dual_hodge3(const CutCellMesh<3> &ccm);
}// namespace mandoline::operators
