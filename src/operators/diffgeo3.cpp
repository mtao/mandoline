#include "mandoline/operators/diffgeo3.hpp"
#include "mandoline/operators/boundary3.hpp"
#include "mandoline/operators/volume3.hpp"

namespace mandoline::operators {

// primal-2 form -h> dual-1 -d> dual-0
Eigen::SparseMatrix<double> divergence(const CutCellMesh<3> &ccm) {
    auto B = boundary(ccm,false);
    auto h2 = dual_hodge2(ccm);
    return B * h2.asDiagonal();
}

// primal-3 -d> primal-2 form -h> dual-1 -d> dual-0
Eigen::SparseMatrix<double> laplacian(const CutCellMesh<3> &ccm) {
    auto B = boundary(ccm,false);
    auto h2 = dual_hodge2(ccm);
    return B * h2.asDiagonal() * B.transpose();
}
}
