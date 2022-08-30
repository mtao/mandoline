#pragma once
#include <balsa/eigen/types.hpp>

namespace mandoline::construction::tools {
std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i> remesh_self_intersections(const balsa::eigen::ColVecs3d &V, const balsa::eigen::ColVecs3i &F);
// TODO: pass in the right relative indices
// std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i, balsa::eigen::VecXi,balsa::eigen::VecXi> remesh_self_intersections_with_indices(const balsa::eigen::ColVecs3d &V, const balsa::eigen::ColVecs3i &F);
}
