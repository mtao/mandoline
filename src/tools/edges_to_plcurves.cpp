#include <balsa/eigen/types.hpp>
#include "mandoline/tools/edges_to_plcurves.hpp"
#include "mandoline/construction/face_collapser.hpp"
#include <balsa/eigen/stl2eigen.hpp>
#include <mtao/geometry/mesh/edges_to_plcurves.hpp>
#include <mtao/geometry/volume.hpp>
#include <mtao/geometry/mesh/unique_simplices.hpp>
#include <spdlog/spdlog.h>


namespace mandoline::tools {
std::vector<std::tuple<std::vector<int>, bool>> edge_to_plcurves(
  const balsa::eigen::ColVecs2d &V,
  const balsa::eigen::ColVecs2i &E,
  bool closed_only) {
    construction::FaceCollapser fc(E);
    fc.bake(V, false);// ignore finding faces in faces
    std::vector<std::tuple<std::vector<int>, bool>> ret;
    for (auto &&[fidx, edges] : fc.face_edges()) {

        // if the edge is open then i need to do this to remove pairs of edges
        // this preserves input orientation when possible :)
        balsa::eigen::ColVecs2i E = mtao::geometry::mesh::unique_simplices(balsa::eigen::stl2eigen(edges));
        auto r = mtao::geometry::mesh::edge_to_plcurves(E, closed_only);

        //std::cout << "Fidx: " << fidx << std::endl;
        double ovol = mtao::geometry::brep_volume(V,E);
        r.erase(std::remove_if(r.begin(),r.end(),[&](auto&& pr){
                    auto&& [vec, closed] = pr;
                    //std::copy(vec.begin(),vec.end(),std::ostream_iterator<int>(std::cout,","));
                    if(closed) {
                    double vol = mtao::geometry::curve_volume(V, vec);
                    //std::cout << " vols: " << vol << "/" << ovol;
                    //std::cout << "(" << (vol * ovol < 1e-10) << ")";
                    //std::cout << std::endl;
                    return vol * ovol < -1e-10;

                    } else {
                    //std::cout << std::endl;
                    return false;
                    }
                    }), r.end());

        ret.insert(ret.end(), r.begin(), r.end());
    }
    return ret;
}

}// namespace mandoline::tools
