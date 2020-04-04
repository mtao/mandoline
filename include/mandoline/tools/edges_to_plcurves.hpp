#pragma once
#include <mtao/types.hpp>

namespace mandoline::tools {
    // closed_only makes open curves wrap around themselves
    // two_sided returns edges on both sides of boundary curves
std::vector<std::tuple<std::vector<int>, bool>> edge_to_plcurves(
    const mtao::ColVecs2d& V,
    const mtao::ColVecs2i& E, bool closed_only = true);

}
