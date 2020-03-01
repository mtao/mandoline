#pragma once

#include <mtao/types.hpp>
#include <Eigen/Sparse>
#include "mandoline/cutedge.hpp"

namespace mandoline {
struct InterpolatedEdge {

    // the interpolated coordinates
    mtao::Vec2d ts;
    int parent_fid = -1;

    std::map<std::array<int, 2>, double> sparse_matrix_entries(const std::array<int,2> &indices, const mtao::ColVecs2i &E) const;
    template <int D>
    std::map<std::array<int, 2>, double> sparse_matrix_entries(const CutEdge<D> &edge, const mtao::ColVecs2i &E) const {
        assert(edge.is_mesh_edge());
        return sparse_matrix_entries(*edge.indices.begin(), E);
    }
    std::map<std::array<int, 2>, double> sparse_matrix_entries(const CutMeshEdge<2> &edge, const mtao::ColVecs2i &E) const {
        return sparse_matrix_entries(edge.indices, E);
    }
    double volume() const;
};
}// namespace mandoline
