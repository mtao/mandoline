#pragma once
#include <mtao/types.hpp>

namespace mandoline::construction::tools {
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> remesh_self_intersections(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F);
// TODO: pass in the right relative indices
// std::tuple<mtao::ColVecs3d, mtao::ColVecs3i, mtao::VecXi,mtao::VecXi> remesh_self_intersections_with_indices(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F);
}
