#pragma once
#include "mandoline/mesh2.hpp"


namespace mandoline::operators {
//mesh vertex -> cut vertex
Eigen::SparseMatrix<double> barycentric_matrix(const CutCellMesh<2> &ccm);
//mesh face -> cut face
Eigen::SparseMatrix<double> face_barycentric_volume_matrix(const CutCellMesh<2> &ccm);
//grid vertex -> cut vertex
Eigen::SparseMatrix<double> trilinear_matrix(const CutCellMesh<2> &ccm);
//grid face -> cut face
Eigen::SparseMatrix<double> face_grid_volume_matrix(const CutCellMesh<2> &ccm);
}// namespace mandoline::operators
