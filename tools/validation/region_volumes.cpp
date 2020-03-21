#define __gmp_const const
#include "cutmesh_validation.hpp"
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Gmpq.h>
#include <igl/copyleft/cgal/extract_cells.h>
#include <mtao/geometry/bounding_box_mesh.hpp>
#include <set>
using namespace mandoline;
bool equal_region_check(const mandoline::CutCellMesh<3> &ccm) {
    auto C = input_mesh_regions(ccm);
    auto [a, b] = region_counts(ccm, C);
    return a == b;
}

std::array<int, 2> region_counts(const mandoline::CutCellMesh<3> &ccm) {
    auto C = input_mesh_regions(ccm);
    return region_counts(ccm, C);
}

std::array<int, 2> region_counts(const mandoline::CutCellMesh<3> &ccm, const mtao::ColVecs2i &C) {
    auto ccm_regions = ccm.regions();
    //igl's C will give us regions, including the infinite region.
    //|regions - {infinite}| = maxCoeff because of an implicit "-1" for region 0
    std::set<int> cr;
    std::copy(ccm_regions.begin(), ccm_regions.end(), std::inserter(cr, cr.end()));
    return { { int(cr.size()), C.maxCoeff() } };
}

Eigen::VectorXd brep_region_volumes(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F) {
    auto C = regions(V, F);
    return brep_region_volumes(V, F, C);
}
Eigen::VectorXd brep_region_volumes(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const mtao::ColVecs2i &C) {

    int cells = C.maxCoeff();
    Eigen::VectorXd vol = Eigen::VectorXd::Zero(cells + 1);
    for (int i = 0; i < F.cols(); ++i) {
        Eigen::Matrix3d A;
        auto f = F.col(i);
        auto c = C.col(i);
        for (int j = 0; j < 3; ++j) {
            A.col(j) = V.col(f(j));
        }
        double v = A.determinant();
        vol(c(0)) -= v;
        vol(c(1)) += v;
    }
    vol /= 6;
    return vol;
}
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> orig_with_bbox(const mandoline::CutCellMesh<3> &ccm) {
    auto &&V = ccm.origV();
    auto &&F = ccm.origF();

    auto [gV, gF] = mtao::geometry::bounding_box_mesh(ccm.bbox());
    mtao::ColVecs3d RV = mtao::eigen::hstack(V, gV);
    mtao::ColVecs3i RF = mtao::eigen::hstack(F, gF.array() + V.cols());
    return { RV, RF };
}
mtao::VecXd brep_region_volumes(const mandoline::CutCellMesh<3> &ccm) {
    auto [V, F] = orig_with_bbox(ccm);
    auto C = input_mesh_regions(ccm);
    return brep_region_volumes(V, F, C).bottomRows(C.size() - 1);
}

std::map<int, double> region_volumes(const mandoline::CutCellMesh<3> &ccm) {
    std::map<int, int> unindexer;
    auto R = ccm.regions();
    auto V = ccm.cell_volumes();
    for (int i = 0; i < R.size(); ++i) {
        if (unindexer.find(R[i]) == unindexer.end()) {
            int s = unindexer.size();
            unindexer[R[i]] = s;
        }
    }

    mtao::VecXd vol(unindexer.size());
    vol.setZero();
    assert(R.size() == V.size());
    for (int i = 0; i < ccm.cell_size(); ++i) {
        vol(unindexer[R[i]]) += V(i);
    }
    std::map<int, double> vols;
    for (auto &&[a, b] : unindexer) {
        vols[b] = vol(a);
    }
    return vols;
}


mtao::ColVecs2i input_mesh_regions(const mandoline::CutCellMesh<3> &ccm) {
    Eigen::MatrixXi C;
    auto V = ccm.origV();
    auto F = ccm.origF();
    return regions(V, F);
}

mtao::ColVecs2i regions(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F) {
    Eigen::MatrixXi C;
    Eigen::MatrixXd VV = V.transpose();
    Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, Eigen::Dynamic, Eigen::Dynamic> VVV = VV.cast<CGAL::Lazy_exact_nt<CGAL::Gmpq>>();
    Eigen::MatrixXi FF = F.transpose();
    igl::copyleft::cgal::extract_cells(VV, FF, C);
    //igl::copyleft::cgal::extract_cells(VV.cast<CGAL::Lazy_exact_nt<CGAL::Gmpq>>().eval(),FF,C);
    return C.transpose();
}
