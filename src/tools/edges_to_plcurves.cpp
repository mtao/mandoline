#include <mtao/types.hpp>
#include "mandoline/tools/edges_to_plcurves.hpp"
#include "mandoline/construction/face_collapser.hpp"
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/geometry/mesh/edges_to_plcurves.hpp>
#include <mtao/geometry/volume.hpp>
#include <spdlog/spdlog.h>


namespace mandoline::tools {
std::vector<std::tuple<std::vector<int>, bool>> edge_to_plcurves(
  const mtao::ColVecs2d &V,
  const mtao::ColVecs2i &E,
  bool closed_only) {
    construction::FaceCollapser fc(E);
    fc.bake(V, false);// ignore finding faces in faces
    std::vector<std::tuple<std::vector<int>, bool>> ret;
    for (auto &&[fidx, edges] : fc.face_edges()) {
        auto r = mtao::geometry::mesh::edge_to_plcurves(mtao::eigen::stl2eigen(edges), closed_only);

        mtao::ColVecs2i E = mtao::eigen::stl2eigen(edges);
        spdlog::info("Num verts: {}, E range [{},{}]", V.cols(), E.minCoeff(),E.maxCoeff());
        double ovol = mtao::geometry::brep_volume(V,E);
        r.erase(std::remove_if(r.begin(),r.end(),[&](auto&& pr){
                    auto&& [vec, closed] = pr;
                    if(closed) {
                    double vol = mtao::geometry::curve_volume(V, vec);
                    return vol * ovol < 1e-10;

                    } else {
                    return false;
                    }
                    }), r.end());

        ret.insert(ret.end(), r.begin(), r.end());
    }
    return ret;
}

}// namespace mandoline::tools
