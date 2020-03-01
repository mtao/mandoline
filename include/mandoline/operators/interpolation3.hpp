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
Eigen::SparseMatrix<double> face_grid_volume_matrix(const CutCellMesh<3> &ccm);

}// namespace mandoline::operators
