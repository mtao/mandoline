#include "mandoline/barycentric_triangle_face.hpp"
#include <iostream>
#include <mtao/geometry/volume.h>
namespace mandoline {
std::map<std::array<int, 2>, double> BarycentricTriangleFace::sparse_matrix_entries(const std::vector<int> &indices, const balsa::eigen::ColVecs3i &F) const {
    auto f = F.col(parent_fid);
    std::map<std::array<int, 2>, double> ret;
    assert(indices.size() == barys.cols());
    for (int i = 0; i < indices.size(); ++i) {
        auto B = barys.col(i);
        auto row = indices[i];
        for (int j = 0; j < 3; ++j) {
            if (B(j) != 0) {
                //ret.emplace_back(row,f(j),B(j));
                ret[std::array<int, 2>{ { row, f(j) } }] = B(j);
            }
        }
    }
    return ret;
}
double BarycentricTriangleFace::volume() const {
    return mtao::geometry::curve_volume(barys.topRows<2>());
}
}// namespace mandoline
