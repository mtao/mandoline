#include <mtao/types.hpp>
#include "mandoline/tools/edges_to_plcurves.hpp"
#include "mandoline/construction/face_collapser.hpp"
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/geometry/mesh/edges_to_plcurves.hpp>


namespace mandoline::tools {
std::vector<std::tuple<std::vector<int>, bool>> edge_to_plcurves(
  const mtao::ColVecs2d &V,
  const mtao::ColVecs2i &E,
  bool closed_only) {
    construction::FaceCollapser fc(E);
    fc.bake(V, false);// ignore finding faces in faces
    std::vector<std::tuple<std::vector<int>, bool>> ret;
    for (auto &&[fidx, edges] : fc.face_edges()) {
        std::cout << mtao::eigen::stl2eigen(edges) << std::endl;
        auto r = mtao::geometry::mesh::edge_to_plcurves(mtao::eigen::stl2eigen(edges), closed_only);
        ret.insert(ret.end(), r.begin(), r.end());
    }
    return ret;
}

}// namespace mandoline::tools
