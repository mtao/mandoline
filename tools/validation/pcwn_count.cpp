#define __gmp_const const
#include <igl/copyleft/cgal/piecewise_constant_winding_number.h>
#include "cutmesh_validation.hpp"
using namespace mandoline;
std::array<int,2> pcwn_count(const CutCellMesh<3>& ccm) {
    auto V = ccm.vertices();
    Eigen::MatrixXd IV = V.transpose();
    int is_pcwn = 0;
    int is_not_pcwn = 0;
    for(auto&& cid: mtao::iterator::range(ccm.cells().size())) {
        auto [V_,F] = ccm.triangulated_cell(true,true);
        mtao::RowVecs3i FF = F.transpose();
        double wn = 0;
        if(V_.size() == 0) {
            wn = igl::copyleft::cgal::piecewise_constant_winding_number(IV,FF);
        } else {
            mtao::RowVecs3d VV = V_.transpose();
            wn = igl::copyleft::cgal::piecewise_constant_winding_number(VV,FF);
        }

        if(wn) {
            is_pcwn++;
        } else {
            mtao::logging::warn() << "Not pcwn cell: " << cid;
            is_not_pcwn++;
        }
    }
    return {{is_pcwn,is_not_pcwn}};
}
bool pcwn_check(const CutCellMesh<3>& ccm) {
    return std::get<1>(pcwn_count(ccm)) == 0;
}
