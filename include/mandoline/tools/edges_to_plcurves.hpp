#pragma once
#include <mtao/types.hpp>

namespace mandoline::tools {
std::vector<std::tuple<std::vector<int>, bool>> edge_to_plcurves(
    const mtao::ColVecs2d& V,
    const mtao::ColVecs2i& E, bool closed_only = true);

}
