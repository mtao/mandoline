#pragma once
#include <mtao/types.hpp>

namespace mandoline::construction {
// the internal core of mandoline assumes that the faces are unique, so this removes self intersections and makes the entries unique.
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> preprocess_mesh(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F);
}
