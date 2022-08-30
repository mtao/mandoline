#include "mandoline/construction/tools/remesh_self_intersections.hpp"
#ifdef MANDOLINE_HANDLE_SELF_INTERSECTIONS
#define __gmp_const const
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/remove_unreferenced.h>
#endif

#include <mtao/logging/profiler.hpp>

namespace mandoline::construction::tools {
std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i> remesh_self_intersections(const balsa::eigen::ColVecs3d &V, const balsa::eigen::ColVecs3i &F) {
    auto t = mtao::logging::profiler("Remeshing self intersections");
#ifdef MANDOLINE_HANDLE_SELF_INTERSECTIONS
    Eigen::MatrixXd VV, SV;
    Eigen::MatrixXi FF, SF, IF;
    Eigen::VectorXi J, IM, UIM;

    igl::copyleft::cgal::remesh_self_intersections(V.transpose(), F.transpose(), { false, false, true }, VV, FF, IF, J, IM);
    std::for_each(FF.data(), FF.data() + FF.size(), [&IM](int &a) { a = IM(a); });
    igl::remove_unreferenced(VV, FF, SV, SF, UIM);
    return std::make_tuple(balsa::eigen::ColVecs3d(SV.transpose()), balsa::eigen::ColVecs3i(SF.transpose()));
#else
    mtao::logging::error() << "mandoline::construction::remesh_self_intersections was called even though it was built without the ability to handle such things! Try building with ``-DHANDLE_SELF_INTERSECTIONS";
#endif//MANDOLINE_HANDLE_SELF_INTERSECTIONS
    return { V, F };
}
}// namespace mandoline::construction
