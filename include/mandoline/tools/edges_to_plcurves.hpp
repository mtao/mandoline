#pragma once
#include <balsa/eigen/types.hpp>

namespace mandoline::tools {
    // closed_only makes open curves wrap around themselves
    // two_sided returns edges on both sides of boundary curves
std::vector<std::tuple<std::vector<int>, bool>> edge_to_plcurves(
    const balsa::eigen::ColVecs2d& V,
    const balsa::eigen::ColVecs2i& E, bool closed_only = true);

}
