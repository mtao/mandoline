#pragma once


#include <mandoline/mesh3.hpp>



std::map<int,int> face_regions(const mandoline::CutCellMesh<3>& ccm);
mandoline::CutCellMesh<3> read(const std::string& filename);
mtao::VecXd flux(const mandoline::CutCellMesh<3>& ccm);
mtao::VecXd divergence(const mandoline::CutCellMesh<3>& ccm);
Eigen::SparseMatrix<double> boundary(const mandoline::CutCellMesh<3>& ccm);
Eigen::SparseMatrix<double> laplacian(const mandoline::CutCellMesh<3>& ccm);
