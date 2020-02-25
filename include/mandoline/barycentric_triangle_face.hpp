#pragma once

#include <mtao/types.hpp>
#include <Eigen/Sparse>
#include "mandoline/cutface.hpp"

namespace mandoline {
struct BarycentricTriangleFace {

    mtao::ColVecs3d barys;
    int parent_fid = -1;

    std::map<std::array<int, 2>, double> sparse_entries(const std::vector<int> &indices, const mtao::ColVecs3i &F) const;
    std::map<std::array<int, 2>, double> sparse_entries(const CutFace<3> &face, const mtao::ColVecs3i &F) const {
        assert(face.is_mesh_face());
        assert(face.indices.size() == 1);
        return sparse_entries(*face.indices.begin(), F);
    }
    std::map<std::array<int, 2>, double> sparse_entries(const CutMeshFace<3> &face, const mtao::ColVecs3i &F) const {
        return sparse_entries(face.indices, F);
    }
    double volume() const;
};
}// namespace mandoline
