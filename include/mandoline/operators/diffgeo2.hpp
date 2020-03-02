#pragma once
#include "mandoline/mesh2.hpp"

namespace mandoline::operators {

// primal-1 form -h> dual-1 -d> dual-0
Eigen::SparseMatrix<double> divergence(const CutCellMesh<2> &ccm);

// primal-2 -d> primal-1 form -h> dual-1 -d> dual-0
Eigen::SparseMatrix<double> laplacian(const CutCellMesh<2> &ccm);

// mesh cut-vertex -> {mesh cut-edge,mesh cut-edge}
std::map<int,std::array<int,2>>  surface_adjacency(const CutCellMesh<2>& ccm);

mtao::VecXd surface_dual_lengths(const CutCellMesh<2>& ccm);
Eigen::SparseMatrix<double> surface_divergence(const CutCellMesh<2>& ccm);
Eigen::SparseMatrix<double> surface_laplacian(const CutCellMesh<2>& ccm);

mtao::VecXd surface_dual_lengths(const CutCellMesh<2>& ccm, const std::map<int,std::array<int,2>>& adj_struct);
Eigen::SparseMatrix<double> surface_boundary(const CutCellMesh<2>& ccm, const std::map<int,std::array<int,2>>& adj_struct);
Eigen::SparseMatrix<double> surface_divergence(const CutCellMesh<2>& ccm, const std::map<int,std::array<int,2>>& adj_struct);
Eigen::SparseMatrix<double> surface_laplacian(const CutCellMesh<2>& ccm, const std::map<int,std::array<int,2>>& adj_struct);
}
