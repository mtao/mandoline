#include "mandoline/operators/diffgeo2.hpp"
#include "mandoline/operators/boundary2.hpp"
#include "mandoline/operators/volume2.hpp"
#include <spdlog/spdlog.h>

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
                adj[edge.indices[0]]= std::array<int, 2>{ { eidx, -1 } };
            } else {
                //it->second[1] = eidx;
                adj[edge.indices[0]][1] = eidx;
            }
            if (auto it = adj.find(edge.indices[1]); it == adj.end()) {
                adj[edge.indices[1]] = std::array<int, 2>{ { eidx, -1 } };
            } else {
                //it->second[0] = eidx;
                adj[edge.indices[1]][1] = eidx;
            }
        }
    }

    // create an orientation by running around open edges and then random edges
    std::set<int> unseen_edges;
    std::set<int> open_edges;
    for(auto&& [k,v]: adj) {
        unseen_edges.emplace(k);
        if(v[0] < 0 || v[1] < 0) {
            open_edges.emplace(k);
        }
    }
    auto remove = [&](int idx) {

        if(auto it = open_edges.find(idx); it != open_edges.end()) {
            open_edges.erase(it);
        }
        if(auto it = unseen_edges.find(idx); it != unseen_edges.end()) {
            unseen_edges.erase(it);
        }
    };
    auto run = [&](const int start_vertex, const int start_edge) {

        int vertex = start_vertex;
        int edge = start_edge;
        do {
            auto&& e = ccm.cut_edge(start_edge).indices;
            // step vertex
            if(e[0] == vertex) {
                vertex = e[1];
            } else {
                vertex = e[0];
            }
            remove(vertex);

            // step edge and swap it into the right order
            auto& de = adj[vertex];
            if(de[0] != edge) {
                std::swap(de[0],de[1]);
            }
            edge = de[1];
            if(edge < 0) {
                break;
            }


        } while(vertex != start_vertex);

    };
    while(!open_edges.empty()) {

        int start_idx = *open_edges.begin();
        remove(start_idx);
        int start_edge;
        auto&& e = adj[start_idx];
        if(e[0] == -1) {
            start_edge = e[1];
        } else {
            start_edge = e[0];
        }
        run(start_idx,start_edge);
    }
    while(!unseen_edges.empty()) {
        int start_idx = *unseen_edges.begin();
        remove(start_idx);
        auto&& e = adj[start_idx];
        int start_edge = e[0];

        run(start_idx,start_edge);

    }
    return adj;
}
balsa::eigen::VecXd surface_dual_lengths(const CutCellMesh<2> &ccm) {
    return surface_dual_lengths(ccm, surface_adjacency(ccm));
}
Eigen::SparseMatrix<double> surface_boundary(const CutCellMesh<2> &ccm) {
    return surface_boundary(ccm, surface_adjacency(ccm));
}
Eigen::SparseMatrix<double> surface_divergence(const CutCellMesh<2> &ccm) {
    return surface_divergence(ccm, surface_adjacency(ccm));
}
Eigen::SparseMatrix<double> surface_laplacian(const CutCellMesh<2> &ccm) {
    return surface_laplacian(ccm, surface_adjacency(ccm));
}

balsa::eigen::VecXd surface_dual_lengths(const CutCellMesh<2> &ccm, const std::map<int, std::array<int, 2>> &adj_struct) {

    auto vols = edge_lengths(ccm);
    //balsa::eigen::VecXd el(ccm.cut_edges().size());
    balsa::eigen::VecXd el(ccm.num_vertices());
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
            triplets.emplace_back(p, dual_edge_index, 1);
        }
    }
    Eigen::SparseMatrix<double> R(ccm.cut_edges().size(), ccm.vertex_size());
    R.setFromTriplets(triplets.begin(), triplets.end());
    return R;
}
Eigen::SparseMatrix<double> surface_divergence(const CutCellMesh<2> &ccm, const std::map<int, std::array<int, 2>> &adj_struct) {
    balsa::eigen::VecXd h = surface_dual_lengths(ccm, adj_struct);
    auto B = surface_boundary(ccm, adj_struct);
    return B * h.asDiagonal();
}
Eigen::SparseMatrix<double> surface_laplacian(const CutCellMesh<2> &ccm, const std::map<int, std::array<int, 2>> &adj_struct) {

    balsa::eigen::VecXd h = surface_dual_lengths(ccm, adj_struct);
    auto B = surface_boundary(ccm, adj_struct);
    return B * h.asDiagonal() * B.transpose();
}
}// namespace mandoline::operators
