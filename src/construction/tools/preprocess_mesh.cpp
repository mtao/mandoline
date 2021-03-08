#include "mandoline/construction/tools/remesh_self_intersections.hpp"
#include <mtao/geometry/mesh/unique_simplices.hpp>
#include <mtao/geometry/prune_vertices.hpp>

namespace mandoline::construction::tools {

std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> preprocess_mesh(const mtao::ColVecs3d &V_, const mtao::ColVecs3i &F_) {
    mtao::ColVecs3d V;
    mtao::ColVecs3i F;
    std::tie(V, F) = mtao::geometry::prune(V_, F_, 0);
    std::tie(V, F) = remesh_self_intersections(V, F);

    F = mtao::geometry::mesh::unique_simplices(F);
    return { V, F };
}
}// namespace mandoline::construction
