#pragma once

#include <array>
#include <tuple>
#include <mtao/geometry/cyclic_order.hpp>
#include <map>
#include <set>
#include <vector>
#include <spdlog/spdlog.h>
#include <mtao/data_structures/disjoint_set.hpp>
#include <mtao/geometry/mesh/triangle_fan.hpp>
#include <mtao/geometry/trigonometry.hpp>
#include "mandoline/cutface.hpp"
#include <mtao/logging/logger.hpp>
#include <mtao/logging/profiler.hpp>
#include <Eigen/Eigenvalues>


namespace mandoline::construction {
struct CellCollapser {
  public:
    using Edge = std::array<int, 2>;
    using HalfFace = std::tuple<int, Edge>;
    using coord_type = std::array<int, 3>;
    //CellCollapser(const mtao::map<int,std::set<std::vector<int>>>& faces);
    CellCollapser(const mtao::map<int, CutFace<3>> &faces);
    std::tuple<std::set<std::vector<int>>, std::set<Edge>> clean_unmanifold(const std::vector<int> &);
    // for each edge, list off all of the ordered halffaces that use it 
    mtao::map<Edge, std::set<const HalfFace *>> collect_halffaces() const;

    int dual_cell(const HalfFace &hf) const;
    int cell(const HalfFace &hf) const;

    void fill_cell_boundaries();
    // merge_open_faces makes both sides of a triangle part of the same cell if it has an unshared edge
    template<typename Derived>
    void merge(const Eigen::MatrixBase<Derived> &V, bool merge_open_faces = true);
    template<typename Derived>
    void merge_around_edge(const Eigen::MatrixBase<Derived> &V, const Edge& e, const std::set<const HalfFace *>& halffaces);
    template<typename Derived>
    void bake(const Eigen::MatrixBase<Derived> &V);

    std::vector<std::set<int>> cell_faces() const;
    const auto &faces() const { return m_faces; }
    const auto &face(int idx) const { return m_faces.at(idx); }

    const auto &cell_boundaries() const { return m_cell_boundaries; }
    void remove_boundary_cells();
    void remove_grid_boundary_cells(const std::array<int, 3> &shape);
    void remove_boundary_cells_by_volume(const std::map<int, double> &face_vols);
    void remove_boundary_cells_from_vertices(const std::set<int> &boundary_vertices);
    void remove_boundary_cells_from_faces(const std::set<int> &boundary_faces);
    std::set<int> folded_faces() const;//faces that have their duals as well
    //private:
    // stores a halfface (edge + face index) and maps it to a cell (the sign maintains the chirality wrt the original face)
    mtao::map<HalfFace, std::tuple<int, bool>> m_halfface_to_cell;
    mtao::map<int, CutFace<3>> m_faces;
    std::map<int, std::set<Edge>> flap_edges;
    mtao::map<int, std::set<coord_type>> face_cell_possibilities;
    std::vector<mtao::map<int, bool>> m_cell_boundaries;
    std::set<Edge> invalid_edges;
    mtao::data_structures::DisjointSet<int> cell_ds;

};

template<typename Derived>
void CellCollapser::bake(const Eigen::MatrixBase<Derived> &V) {
    merge(V);
}

template<typename Derived>
void CellCollapser::merge_around_edge(const Eigen::MatrixBase<Derived> &V, const Edge& e, const std::set<const HalfFace *>& halffaces) {
    //auto t3 = mtao::logging::profiler("cell collapser edge fan processessing",false,"profiler");
    if (halffaces.empty()) {
        return;
    }
    auto [a, b] = e;
    auto va = V.col(a);
    auto vb = V.col(b);
    auto vba = va - vb;

    int maxcoeff;
    // if the edge is really small then we use eigen analysis of the face to compute a normal
    {
        if (vba.norm() < 1e-8) {
            mtao::Mat3d A = mtao::Mat3d::Zero();

            for (auto &&hfp : halffaces) {
                auto [fidx, e] = *hfp;

                bool sign = std::get<1>(m_halfface_to_cell.at(*hfp));
                if (sign) {
                    auto &&F = m_faces.at(fidx);
                    auto &&N = F.N;
                    A += N * N.transpose();
                }
            }
            Eigen::SelfAdjointEigenSolver<mtao::Mat3d> eigensolver(A);
            if (eigensolver.info() != Eigen::Success) {
                vba.cwiseAbs().maxCoeff(&maxcoeff);
            } else {
                eigensolver.eigenvectors().col(0).cwiseAbs().maxCoeff(&maxcoeff);
            }
        } else {
            vba.cwiseAbs().maxCoeff(&maxcoeff);
        }
    }
    const bool flip_axes = vba(maxcoeff) < 0;

    std::vector<int> faces(halffaces.size());
    mtao::ColVecs2d A(2,halffaces.size());

    for (auto && [idx,fidx,hfp] : mtao::iterator::enumerate(faces,halffaces)) {
        auto [fidx_, e] = *hfp;

        fidx = fidx_;


        auto &&F = m_faces.at(fidx);
        const bool sign = std::get<1>(m_halfface_to_cell.at(*hfp));
        mtao::Vec3d N = (sign ? 1 : -1) * F.N.normalized();

        auto a = A.col(idx);
        if(flip_axes) {
            a(0) = N((maxcoeff+1)%3);
            a(1) = N((maxcoeff+2)%3);
        } else {
            a(1) = N((maxcoeff+1)%3);
            a(0) = N((maxcoeff+2)%3);
        }

    }
    auto order = mtao::geometry::cyclic_order(A);
    std::vector<int> ordered_faces(halffaces.size());
    // dereference the indices of indices
    std::transform(order.begin(), order.end(), ordered_faces.begin(), [&](int idx) { return faces[idx]; });
    auto it = ordered_faces.begin();
    auto it1 = it;
    it1++;
    for (; it != ordered_faces.end(); ++it, ++it1) {
        if (it1 == ordered_faces.end()) {
            it1 = ordered_faces.begin();
        }
        int idx = cell({*it,e});
        int idx1 = dual_cell({*it1,e});
        cell_ds.join(idx,idx1);
    }
}
template<typename Derived>
void CellCollapser::merge(const Eigen::MatrixBase<Derived> &V, bool merge_open_faces) {
    auto t2 = mtao::logging::profiler("cell collapser merge", false, "profiler");
    auto&& hfs = collect_halffaces();
    for (auto [e, halffaces] : hfs) {
        if (e[0] < e[1]) { // only need to process each edge once
            if(merge_open_faces || halffaces.size() > 1) {
                merge_around_edge(V,e,halffaces);
            }
        }
    }
    fill_cell_boundaries();
}
}// namespace mandoline::construction
