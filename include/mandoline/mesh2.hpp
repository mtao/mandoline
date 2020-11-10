#pragma once
#include "mandoline/mesh.hpp"
#include "mandoline/cutface.hpp"
#include "mandoline/exterior_grid.hpp"
namespace mandoline {
template<>
struct CutCellMesh<2> : public CutCellMeshBase<2, CutCellMesh<2>> {
  public:
    using Base = CutCellMeshBase<2, CutCellMesh<2>>;
    using Base::Base;
    using ColVecs = typename Base::ColVecs;
    using VecX = typename Base::VecX;

    int num_cells() const;
    int num_cutcells() const;
    int num_faces() const;
    int num_cutfaces() const;
    int num_edges() const;
    int num_cutedges() const;
    int cell_index(const VecCRef &) const;
    coord_type grid_cell_coord(const VecCRef &p) const { return std::get<0>(StaggeredGrid::coord(p)); }
    int grid_cell_index(const VecCRef &p) const { return cell_grid().index(grid_cell_coord(p)); }
    int grid_cell_index(const coord_type &c) const { return cell_grid().index(c); }

    mtao::ColVecs2i dual_edges() const;
    Edges edges() const;

    //ColVecs dual_vertices() const;
    mtao::ColVectors<int, 3> faces() const;
    ColVecs centroids() const;
    std::set<std::vector<int>> cell(int index) const;

    int nearest_edge_index(const VecCRef &) const;
    bool in_cell(const VecCRef &, int idx) const;
    bool in_cell(const ColVecs &V, const VecCRef &, int idx) const;
    VecX volumes() const;
    VecX dual_edge_volumes() const;
    mtao::VecXi regions() const;
    const std::map<int, std::map<int, bool>> &face_boundary_map() const { return m_face_boundary_map; }
    Eigen::SparseMatrix<double> boundary(bool include_domain_boundary_faces) const;
    const std::vector<CutFace<2>> &cut_faces() const { return m_faces; }
    bool is_cutface_index(int index) const { return index < num_cutfaces(); }
    Edges halfedges_per_edge;//halfedge indices

    std::vector<CutFace<2>> m_faces;
    std::map<int, std::map<int, bool>> m_face_boundary_map;
    mtao::geometry::mesh::HalfEdgeMesh hem;
    ExteriorGrid<2> exterior_grid;
    VecX halfedge_orientations;
    std::map<int, int> halfedge_to_edge_index;
    std::map<int, std::set<int>> cell_grid_ownership;
    //std::map<int, std::set<int>> edge_grid_ownership;
    //std::map<Edge, int> edge_vertices_to_edge_index;

    // should use exterior_grid instead
};
}// namespace mandoline
