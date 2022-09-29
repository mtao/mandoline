#pragma once

#include <balsa/eigen/types.hpp>
#include <Eigen/Sparse>
#include "mandoline/cutface.hpp"

namespace mandoline {
    // BarycentricTriangleFan stores the vertices of a simple polygon
struct BarycentricTriangleFace {

    balsa::eigen::ColVecs3d barys;
    int parent_fid = -1;

    // reads off the vertex info in the format {cut-vertex index, triangle-vertex} = barycentric
    // this format is creates 
    std::map<std::array<int, 2>, double> sparse_matrix_entries(const std::vector<int> &indices, const balsa::eigen::ColVecs3i &F) const;
    // we can just pass a cut-face in. Cut-faces of trianglees only have a single face so we can just extract that one
    std::map<std::array<int, 2>, double> sparse_matrix_entries(const CutFace<3> &face, const balsa::eigen::ColVecs3i &F) const {
        assert(face.is_mesh_face());
        assert(face.indices.size() == 1);
        return sparse_matrix_entries(*face.indices.begin(), F);
    }


    std::map<std::array<int, 2>, double> sparse_matrix_entries(const CutMeshFace<3> &face, const balsa::eigen::ColVecs3i &F) const {
        return sparse_matrix_entries(face.indices, F);
    }
    double volume() const;
};
}// namespace mandoline
