#pragma once

#include <mtao/types.hpp>
#include <Eigen/Sparse>

// TODO: implement this instead of interpolated
namespace mandoline {
    template <int SimplexDim, int ChildPolygonVertexCount = Eigen::Dynamic>
struct SimplexChildPolygon {

    mtao::Matrix<double,SimplexDim,ChildPolygonVertexCount> coordinates;
    int parent_id = -1;

    // reads off the vertex info in the format {cut-vertex index, triangle-vertex} = barycentric
    // this format is creates 
    template <typename PolygonContainer>
    std::vector<Eigen::Triplet<double>> sparse_matrix_entries(const PolygonIndexContainer &indices, const mtao::ColVectors<double,SimplexDim> &S) const;

    template <typename PolygonContainer>
    std::map<std::array<int, 2>, double> sparse_matrix_entries(const CoordMaskedGeometry<D,PolygonContainer> &pg, const Eigen::MatrixBase<Derived> &S) const {
        return sparse_matrix_entries(pg.indices, S);
    }

    // no default impl exists, others need to impl it
    double volume() const;
};
}// namespace mandoline
#include "mandoline/simplex_child_polygon_impl.hpp"
