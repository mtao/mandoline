#pragma once
#include "mandoline/mesh2.hpp"


namespace mandoline::operators {
balsa::eigen::VecXd face_volumes(const CutCellMesh<2> &ccm, bool from_triangulation = false);
balsa::eigen::VecXd edge_lengths(const CutCellMesh<2> &ccm);
balsa::eigen::VecXd dual_edge_lengths(const CutCellMesh<2> &ccm);


// o^k/|a^k| = o^{n-k} / |a^{n-k}|
// H: o^k -> o^{n-k} = |a^{n-k}|/|a^k|

// dual_edge_lengths / face_volumes
// dual 1 -> primal 2
balsa::eigen::VecXd primal_hodge1(const CutCellMesh<2> &ccm);

// face_volumes / dual_edge_lengths
// primal 1 -> dual 1
balsa::eigen::VecXd dual_hodge1(const CutCellMesh<2> &ccm);

// dual_edge_lengths / face_volumes
// dual 0 -> primal 2
balsa::eigen::VecXd primal_hodge2(const CutCellMesh<2> &ccm);

// face_volumes / dual_edge_lengths
// primal 2 -> dual 0
balsa::eigen::VecXd dual_hodge2(const CutCellMesh<2> &ccm);

}// namespace mandoline::operators
