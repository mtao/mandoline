#pragma once
#include "mandoline/mesh3.hpp"

namespace mandoline::operators {

// primal-1 form -h> dual-1 -d> dual-0
Eigen::SparseMatrix<double> divergence(const CutCellMesh<3> &ccm);

// primal-3 -d> primal-1 form -h> dual-1 -d> dual-0
Eigen::SparseMatrix<double> laplacian(const CutCellMesh<3> &ccm);
}
