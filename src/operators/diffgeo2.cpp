#include "mandoline/operators/diffgeo2.hpp"
#include "mandoline/operators/boundary2.hpp"
#include "mandoline/operators/volume2.hpp"

namespace mandoline::operators {

// primal-1 form -h> dual-1 -d> dual-0
Eigen::SparseMatrix<double> divergence(const CutCellMesh<2> &ccm, bool include_domain_boundary) {
    auto B = boundary(ccm, include_domain_boundary);
    auto h1 = dual_hodge1(ccm);
    return B.transpose() * h1.asDiagonal();
}

// primal-2 -d> primal-1 form -h> dual-1 -d> dual-0
Eigen::SparseMatrix<double> laplacian(const CutCellMesh<2> &ccm) {
    auto B = boundary(ccm, false);
    auto h1 = dual_hodge1(ccm);
    return B.transpose() * h1.asDiagonal() * B;
}

std::map<int, std::array<int, 2>> surface_adjacency(const CutCellMesh<2> &ccm) {
    std::map<int, std::array<int, 2>> adj;
    for (auto &&[eidx, edge] : mtao::iterator::enumerate(ccm.cut_edges())) {
        if (edge.is_mesh_edge()) {
            if (auto it = adj.find(edge.indices[0]); it == adj.end()) {
                it->second = std::array<int, 2>{ { -1, eidx } };
            } else {
                it->second[1] = eidx;
            }
            if (auto it = adj.find(edge.indices[1]); it == adj.end()) {
                it->second = std::array<int, 2>{ { eidx, -1 } };
            } else {
                it->second[0] = eidx;
            }
        }
    }
    return adj;
}
mtao::VecXd surface_dual_lengths(const CutCellMesh<2> &ccm) {
    return surface_dual_lengths(ccm, surface_adjacency(ccm));
}
Eigen::SparseMatrix<double> surface_divergence(const CutCellMesh<2> &ccm) {
    return surface_divergence(ccm, surface_adjacency(ccm));
}
Eigen::SparseMatrix<double> surface_laplacian(const CutCellMesh<2> &ccm) {
    return surface_laplacian(ccm, surface_adjacency(ccm));
}

mtao::VecXd surface_dual_lengths(const CutCellMesh<2> &ccm, const std::map<int, std::array<int, 2>> &adj_struct) {

    auto vols = edge_lengths(ccm);
    mtao::VecXd el(ccm.cut_edges().size());
    el.setZero();
    for (auto &&[dual_edge_index, pr] : adj_struct) {
        for (auto &&idx : pr) {
            if (idx >= 0) {
                el(dual_edge_index) += vols(idx);
            }
        }
    }
    el /= 2;
    return el;
}
Eigen::SparseMatrix<double> surface_boundary(const CutCellMesh<2> &ccm, const std::map<int, std::array<int, 2>> &adj_struct) {
    std::vector<Eigen::Triplet<double>> triplets;

    for (auto &&[dual_edge_index, pr] : adj_struct) {
        auto [n, p] = pr;
        if (n >= 0) {
            triplets.emplace_back(n, dual_edge_index, -1);
        }
        if (p >= 0) {
            triplets.emplace_back(p, dual_edge_index, -1);
        }
    }
    Eigen::SparseMatrix<double> R(ccm.cut_edges().size(), ccm.vertex_size());
    R.setFromTriplets(triplets.begin(), triplets.end());
    return R;
}
Eigen::SparseMatrix<double> surface_divergence(const CutCellMesh<2> &ccm, const std::map<int, std::array<int, 2>> &adj_struct) {
    mtao::VecXd h = surface_dual_lengths(ccm, adj_struct);
    auto B = surface_boundary(ccm, adj_struct);
    return B * h.asDiagonal();
}
Eigen::SparseMatrix<double> surface_laplacian(const CutCellMesh<2> &ccm, const std::map<int, std::array<int, 2>> &adj_struct) {

    mtao::VecXd h = surface_dual_lengths(ccm, adj_struct);
    auto B = surface_boundary(ccm, adj_struct);
    return B * h.asDiagonal() * B.transpose();
}
}// namespace mandoline::operators
