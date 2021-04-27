#if !defined(__gmp_const)
#define __gmp_const const
#endif
#ifdef MANDOLINE_ASSUME_CELLS_DONT_INTERSECT
#include <igl/piecewise_constant_winding_number.h>
#else
#include <igl/copyleft/cgal/piecewise_constant_winding_number.h>
#endif
#include <tbb/parallel_reduce.h>

#include "cutmesh_validation.hpp"
using namespace mandoline;
std::array<int, 2> pcwn_count(const CutCellMesh<3>& ccm) {
    auto V = ccm.vertices();

    Eigen::MatrixXd IV = V.transpose();
    return tbb::parallel_reduce(
        tbb::blocked_range<int>(0, ccm.cells().size()),
        std::array<int, 2>{{0, 0}},
        [&](tbb::blocked_range<int> r,
            const std::array<int, 2>& val) -> std::array<int, 2> {
            int is_pcwn = val[0];
            int is_not_pcwn = val[1];
            for (int cid = r.begin(); cid < r.end(); ++cid) {
                auto [V_, F] = ccm.triangulated_cell(true, true);
                mtao::RowVecs3i FF = F.transpose();
                double wn = 0;
#if defined(MANDOLINE_ASSUME_CELLS_DONT_INTERSECT)
                wn = igl::piecewise_constant_winding_number(FF);
#else
                if (V_.size() == 0) {
                    wn = igl::copyleft::cgal::piecewise_constant_winding_number(
                        IV, FF);
                } else {
                    mtao::RowVecs3d VV = V_.transpose();
                    wn = igl::copyleft::cgal::piecewise_constant_winding_number(
                        VV, FF);
                }
#endif

                if (wn) {
                    is_pcwn++;
                } else {
                    mtao::logging::warn() << "Not pcwn cell: " << cid;
                    is_not_pcwn++;
                }
            }
            return {{is_pcwn, is_not_pcwn}};
        },
        [&](const std::array<int, 2>& a,
            const std::array<int, 2>& b) -> std::array<int, 2> {
            std::array<int, 2> r;
            std::transform(a.begin(), a.end(), b.begin(), r.begin(),
                           std::plus<int>());
            return r;
        });
}
bool pcwn_check(const CutCellMesh<3>& ccm) {
    return std::get<1>(pcwn_count(ccm)) == 0;
}
