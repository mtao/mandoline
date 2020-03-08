#pragma once
#include "mandoline/mesh3.hpp"


namespace mandoline::operators {
//mesh vertex -> cut vertex
Eigen::SparseMatrix<double> barycentric_matrix(const CutCellMesh<3> &ccm);
//mesh face -> cut face
Eigen::SparseMatrix<double> face_barycentric_volume_matrix(const CutCellMesh<3> &ccm);
//grid vertex -> cut vertex
Eigen::SparseMatrix<double> trilinear_matrix(const CutCellMesh<3> &ccm);
//grid face -> cut face
// FWIW its worth noting that because cut-faces can lie on axial planes so we optionally can pass in whether or not to include mesh cutfaces.
// include_mesh_cutfaces = true provides a matrix where each grid face is fully utilized.
// include_mesh_cutfaces = false provides a matrix where we partition cut-faces
//  into grid faces and mesh faces 
Eigen::SparseMatrix<double> face_grid_volume_matrix(const CutCellMesh<3> &ccm, bool include_mesh_cutfaces= true);

}// namespace mandoline::operators
