#pragma once
#include <memory>

#include "mandoline/mesh3.hpp"

namespace mandoline::construction {
CutCellMesh<3> from_grid_unnormalized(const mtao::ColVecs3d &V,
                                      const mtao::ColVecs3i &F,
                                      const std::array<int, 3> &cell_shape,
                                      int adaptive_level = 0,
                                      std::optional<double> threshold = 1e-9);
CutCellMesh<3> from_grid(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F,
                         const mtao::geometry::grid::StaggeredGrid3d &grid,
                         int adaptive_level = 0,
                         std::optional<double> threshold = {});
CutCellMesh<3> from_bbox(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F,
                         const Eigen::AlignedBox<double, 3> &bbox,
                         const std::array<int, 3> &cell_shape,
                         int adaptive_level = 0,
                         std::optional<double> threshold = 1e-9);

template <int D>
class CutCellGenerator;
class DeformingGeometryConstructor {
   public:
    DeformingGeometryConstructor(
        const mtao::ColVecs3d &V, const mtao::ColVecs3i &F,
        const mtao::geometry::grid::StaggeredGrid3d &grid,
        int adaptive_level = 0, std::optional<double> threshold = -1);
    DeformingGeometryConstructor();
    DeformingGeometryConstructor(DeformingGeometryConstructor &&o);
    DeformingGeometryConstructor &operator=(DeformingGeometryConstructor &&o);

    ~DeformingGeometryConstructor();
    void update_vertices(const mtao::ColVecs3d &V,
                         const std::optional<double> &threshold = -1);
    void update_topology(const mtao::ColVecs3i &F);
    // updates the mesh and whatnot without topology restrictions
    void set_mesh(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F,
                  const std::optional<double> &threshold = -1);
    // updates the mesh and whatnot without topology restrictions
    void set_vertices(const mtao::ColVecs3d &V,
                      const std::optional<double> &threshold = -1);
    void update_mesh(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F,
                     const std::optional<double> &threshold = -1);
    void update_mesh(const mtao::ColVecs3i &F);
    void update_grid(const mtao::geometry::grid::StaggeredGrid3d &g);
    void set_adaptivity(int res = 0);
    void bake();
    CutCellMesh<3> emit() const;
    void set_external_facet_tossing(bool);
    bool get_external_facet_tossing() const;

    // checks whether the ccg is set up in a valid way
    bool valid() const;

   private:
    std::unique_ptr<CutCellGenerator<3>> _ccg;
    bool _dirty = true;
};
}  // namespace mandoline::construction
