#pragma once
#include <balsa/eigen/types.hpp>

namespace mandoline::construction::tools {
// the internal core of mandoline assumes that the faces are unique, so this removes self intersections and makes the entries unique.
std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i> preprocess_mesh(const balsa::eigen::ColVecs3d &V, const balsa::eigen::ColVecs3i &F);
}
