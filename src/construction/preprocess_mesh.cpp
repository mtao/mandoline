#include "mandoline/construction/remesh_self_intersections.hpp"
#include <mtao/geometry/mesh/unique_simplices.hpp>

namespace mandoline::construction {

std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> preprocess_mesh(const mtao::ColVecs3d &V_, const mtao::ColVecs3i &F_) {
    mtao::ColVecs3d V;
    mtao::ColVecs3i F;
    std::tie(V, F) = remesh_self_intersections(V_, F_);

    F = mtao::geometry::mesh::unique_simplices(F);
    return { V, F };
}
}// namespace mandoline::construction
