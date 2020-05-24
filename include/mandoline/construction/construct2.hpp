#pragma once
#include "mandoline/mesh2.hpp"


namespace mandoline::construction {
CutCellMesh<2> from_grid_unnormalized(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const std::array<int, 2> &cell_shape, std::optional<double> threshold = 1e-9);
CutCellMesh<2> from_grid(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const mtao::geometry::grid::StaggeredGrid2d &grid, std::optional<double> threshold = {});
CutCellMesh<2> from_bbox(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const Eigen::AlignedBox<double, 2> &bbox, const std::array<int, 2> &cell_shape, std::optional<double> threshold = 1e-9);

template<int D>
class CutCellGenerator;
class DeformingGeometryConstructor2 {
  public:
    DeformingGeometryConstructor2(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const mtao::geometry::grid::StaggeredGrid2d &grid, std::optional<double> threshold = -1);
    DeformingGeometryConstructor2(DeformingGeometryConstructor2 &&o);
    DeformingGeometryConstructor2 &operator=(DeformingGeometryConstructor2 &&o);
    DeformingGeometryConstructor2();

    ~DeformingGeometryConstructor2();
    // updates the mesh and whatnot without topology restrictions
    void set_mesh(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const std::optional<double> &threshold = -1);
    // updates the mesh and whatnot without topology restrictions
    void set_vertices(const mtao::ColVecs2d &V, const std::optional<double> &threshold = -1);
    void update_vertices(const mtao::ColVecs2d &V, const std::optional<double> &threshold = -1);
    void update_topology(const mtao::ColVecs2i &F);
    void update_mesh(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const std::optional<double> &threshold = -1);
    void update_mesh(const mtao::ColVecs2i &F);
    void update_grid(const mtao::geometry::grid::StaggeredGrid2d &g);
    void bake();
    CutCellMesh<2> emit() const;

    bool valid() const;

  private:
    std::unique_ptr<CutCellGenerator<2>> _ccg;
    bool _dirty = true;
};
}// namespace mandoline::construction
