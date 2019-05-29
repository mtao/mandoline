#include "mandoline/construction/remesh_self_intersections.hpp"
#define __gmp_const const
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/remove_unreferenced.h>

#include <mtao/logging/profiler.hpp>

namespace mandoline::construction {
    std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> remesh_self_intersections(const mtao::ColVecs3d& V,const mtao::ColVecs3i& F) {
        auto t = mtao::logging::profiler("Remeshing self intersections");
        Eigen::MatrixXd VV,SV;
        Eigen::MatrixXi FF,SF,IF;
        Eigen::VectorXi J,IM,UIM;

        igl::copyleft::cgal::remesh_self_intersections(V.transpose(),F.transpose(),{false,false,true},VV,FF,IF,J,IM);
        std::for_each(FF.data(),FF.data()+FF.size(),[&IM](int & a){a=IM(a);});
        igl::remove_unreferenced(VV,FF,SV,SF,UIM);
        return std::make_tuple(mtao::ColVecs3d(SV.transpose()),mtao::ColVecs3i(SF.transpose()));
    }
}
