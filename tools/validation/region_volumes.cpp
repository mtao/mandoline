#define __gmp_const const
#include "cutmesh_validation.hpp"
#include <igl/copyleft/cgal/extract_cells.h>
#include "mtao/geometry/bounding_box_mesh.hpp"
using namespace mandoline;
bool equal_region_check(const mandoline::CutCellMesh<3>& ccm) {
    auto C = input_mesh_regions(ccm);
    auto [a,b] = region_counts(ccm,C);
}

std::array<int,2> region_counts(const mandoline::CutCellMesh<3>& ccm) {
    auto C = input_mesh_regions(ccm);
    return region_counts(ccm,C);
}

std::array<int,2> region_counts(const mandoline::CutCellMesh<3>& ccm, const mtao::ColVecs2i& C) {
    auto ccm_regions = ccm.regions();
    //igl's C will give us regions, including the infinite region.
    //|regions - {infinite}| = maxCoeff because of an implicit "-1" for region 0
    return {{ccm_regions.size(),C.maxCoeff()}};
}

std::map<int,double> region_volumes(const mandoline::CutCellMesh<3>& ccm) {
    auto C = input_mesh_regions(ccm);

    Eigen::AlignedBox<double,3> bb

    return region_volumes(ccm,C);
}

Eigen::VectorXd brep_region_volumes(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const mtao::ColVecs2i& C) {

    int cells = C.maxCoeff();
    Eigen::VectorXd vol =Eigen::VectorXd::Zero(cells+1);
    for(int i = 0; i < F.cols(); ++i) {
        Eigen::Matrix3d A;
        auto f = F.col(i);
        auto c = C.col(i);
        for(int j = 0; j < 3; ++j) {
            A.col(j) = V.col(f(j));
        }
        double v = A.determinant();
        vol(c(0)) -= v;
        vol(c(1)) += v;


        

    }
    vol /= 6;
    return vol;
}
mtao::VecXd brep_region_volumes(const mandoline::CutCellMesh<3>& ccm) {
    auto C = input_mesh_regions(ccm);
    return brep_region_volumes(ccm,C);
}

std::map<int,double> region_volumes(const mandoline::CutCellMesh<3>& ccm) {
    std::map<int,int> unindexer;
    auto R = ccm.regions();
    auto V = ccm.cell_volumes();

    assert(R.size() == V.size());
    for(int i = 0; i < ccm.num_cells(); ++i) {
    }
}


mtao::ColVecs2i input_mesh_regions(const mandoline::CutCellMesh<3>&) {
    Eigen::MatrixXi C;
    auto V = ccm.origV();
    auto F = ccm.origF();

    igl::copyleft::cgal::extract_cells(V.cast<CGAL::Lazy_exact_nt<CGAL::Gmpq>>().eval(),F,C);
    return C;
}
