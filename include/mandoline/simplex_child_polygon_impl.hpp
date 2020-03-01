#pragma once
#include "mandoline/simplex_child_polygon.hpp"
#include <iostream>
#include <mtao/geometry/volume.h>
namespace mandoline {
    template <int SimplexDim, int ChildPolygonVertexCount>
        template <typename PolygonIndices, typename Derived>
std::vector<Eigen::Triplet<double>> SimplexChildPolygon<SimplexDim, ChildPolygonVertexCount>::sparse_matrix_entries(const PolygonIndices& indices, const Eigen::MatrixBase<Derived> &s) const {
    if constexpr(Derived::ColsAtCompileTime != 1) {
        assert(parent_id < s.cols();
        return sparse_matrix_entries(begin,end,s);
    } else {
    std::map<std::array<int, 2>, double> ret;
    assert(indices.size() == coordinates.cols());

    for (int i = 0; i < indices.size(); ++i) {
        auto B = coordinates.col(i);
        auto row = indices[i];
        for (int j = 0; j < 3; ++j) {
            if (B(j) != 0) {
                //ret.emplace_back(row,f(j),B(j));
                ret[std::array<int, 2>{ { row, s(j) } }] = B(j);
    
        }
    }
    return ret;
    }
}
}// namespace mandoline
