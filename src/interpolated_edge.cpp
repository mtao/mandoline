#include "mandoline/interpolated_edge.hpp"
#include <iostream>
#include <mtao/geometry/volume.h>
namespace mandoline {
std::map<std::array<int, 2>, double> InterpolatedEdge::sparse_matrix_entries(const std::array<int, 2> &indices, const balsa::eigen::ColVecs2i &E) const {
    auto e = E.col(parent_eid);
    std::map<std::array<int, 2>, double> ret;

    for (int i = 0; i < 2; ++i) {
        auto row = indices[i];
        double t = ts(i);
        if (t != 1) {
            ret[std::array<int, 2>{ { row, e[0] } }] = 1 - t;
        }
        if (t != 0) {
            ret[std::array<int, 2>{ { row, e[1] } }] = t;
        }
    }
    return ret;
}
double InterpolatedEdge::volume() const {
    return std::abs(ts(1) - ts(0));
}
}// namespace mandoline
