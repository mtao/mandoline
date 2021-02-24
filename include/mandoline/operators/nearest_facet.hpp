#pragma once
#include "mandoline/mesh2.hpp"
#include "mandoline/mesh3.hpp"
#include "mtao/geometry/grid/staggered_grid.hpp"

namespace igl {
template <typename DerivedV, int DIM>
class AABB;
}
namespace mandoline::operators {

// int nearest_vertex(const CutCellMesh<2>& ccm, Eigen::Ref<const mtao::Vec2d>
// p); int nearest_edge(const CutCellMesh<2>& ccm, Eigen::Ref<const mtao::Vec2d>
// p); int nearest_face(const CutCellMesh<2>& ccm, Eigen::Ref<const mtao::Vec2d>
// p);
//// same as nearest_face
// int nearest_cell(const CutCellMesh<2>& ccm, Eigen::Ref<const mtao::Vec2d> p)
// {
//    return nearest_face(ccm, p);
//}

struct BoundaryFacetProjector3 : public mtao::geometry::grid::StaggeredGrid3d {
    using Base = mtao::geometry::grid::StaggeredGrid3d;
    BoundaryFacetProjector3(const CutCellMesh<3>& ccm);
    ~BoundaryFacetProjector3();
    std::tuple<int, double> nearest_triangle(
        const Eigen::Ref<const mtao::Vec3d> p) const;
    std::tuple<int, double> nearest_grid_face(
        const Eigen::Ref<const mtao::Vec3d> p) const;

    std::vector<std::tuple<int, double>> nearest_triangles(
        const Eigen::Ref<const mtao::ColVecs3d> p) const;
    std::vector<std::tuple<int, double>> nearest_grid_faces(
        const Eigen::Ref<const mtao::ColVecs3d> p) const;
    std::unique_ptr<igl::AABB<mtao::RowVecs3d, 3>> _aabb;
    mtao::RowVecs3d _V;
    mtao::RowVecs3i _F;
};

struct CellParentMaps3 {
    CellParentMaps3(const CutCellMesh<3>& ccm);
    ~CellParentMaps3();
    BoundaryFacetProjector3 _projector;

    std::vector<std::set<int>> grid_contained_cells;
    std::vector<std::set<int>> grid_contained_faces;
    std::vector<std::set<int>> grid_contained_edges;
    std::vector<std::set<int>> grid_contained_vertices;
    std::vector<std::set<int>> cut_cell_coboundary;

    std::vector<std::set<int>> triangle_contained_faces;
};

/// int nearest_vertex(const CutCellMesh<3>& ccm, Eigen::Ref<const mtao::Vec3d>
/// p); int nearest_edge(const CutCellMesh<3>& ccm, Eigen::Ref<const
/// mtao::Vec3d> p); int nearest_face(const CutCellMesh<3>& ccm,
/// Eigen::Ref<const mtao::Vec3d> p); int nearest_cell(const CutCellMesh<3>&
/// ccm, Eigen::Ref<const mtao::Vec3d> p);
///
/// mtao::VecXi nearest_vertices(const CutCellMesh<3>& ccm,
///                             Eigen::Ref<const mtao::ColVecs3d> p);
/// mtao::VecXi nearest_edges(const CutCellMesh<3>& ccm,
///                          Eigen::Ref<const mtao::ColVecs3d> p);
// mtao::VecXi nearest_faces(const CutCellMesh<3>& ccm,
//                          Eigen::Ref<const mtao::ColVecs3d> p);
// mtao::VecXi nearest_faces(const CutCellMesh<3>& ccm,
//                          const CellParentMaps3& parent_maps,
//                          Eigen::Ref<const mtao::ColVecs3d> p);

mtao::VecXi nearest_mesh_cut_faces(const CutCellMesh<3>& ccm,
                                   Eigen::Ref<const mtao::ColVecs3d> p);
mtao::VecXi nearest_mesh_cut_faces(const CutCellMesh<3>& ccm,
                                   const CellParentMaps3& parent_maps,
                                   Eigen::Ref<const mtao::ColVecs3d> p);

mtao::VecXi nearest_cells(const CutCellMesh<3>& ccm,
                          Eigen::Ref<const mtao::ColVecs3d> p);
mtao::VecXi nearest_cells(const CutCellMesh<3>& ccm,
                          const CellParentMaps3& parent_maps,
                          Eigen::Ref<const mtao::ColVecs3d> p);
}  // namespace mandoline::operators
