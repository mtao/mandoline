#pragma once
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mtao/geometry/grid/grid_data.hpp>
#include <set>
#include "cutmesh.pb.h"
#include <mtao/eigen/stl2eigen.hpp>
#include "mandoline/domain_boundary.hpp"
namespace mandoline {
template<int D>
class CutCellMesh;
namespace construction {
    template<int D>
    class CutCellGenerator;
}


template<int D>
class ExteriorGrid : public mtao::geometry::grid::StaggeredGrid<double, D>
  , DomainBoundary {
  public:
    using Base = mtao::geometry::grid::StaggeredGrid<double, D>;
    using GridDatai = mtao::geometry::grid::GridDataD<int, D>;
    using GridDatab = mtao::geometry::grid::GridDataD<bool, D>;
    using coord_type = std::array<int, D>;
    using Edge = std::array<int, 2>;
    using Vec = typename Base::Vec;
    friend class CutCellMesh<D>;
    friend class construction::CutCellGenerator<D>;
    using DomainBoundary::boundary_facet_pairs;


    ExteriorGrid() = default;
    // takes in a cell that indicates true for cells inside the stencil
    ExteriorGrid(const GridDatab &cell_mask);
    ExteriorGrid(const ExteriorGrid &) = default;
    ExteriorGrid(ExteriorGrid &&) = default;
    ExteriorGrid &operator=(const ExteriorGrid &) = default;
    ExteriorGrid &operator=(ExteriorGrid &&) = default;
    mtao::VecXd face_volumes(bool mask_boundary = false) const;
    mtao::VecXd cell_volumes() const;
    int num_faces() const { return DomainBoundary::boundary_facet_size(); }
    int num_cells() const { return m_cell_coords.size(); }
    void make_faces();
    const coord_type &cell_shape() const { return m_cell_indices.shape(); }

    int get_cell_index(const Vec &p) const;


    const GridDatai &cell_indices() const { return m_cell_indices; }

    const std::vector<int> &boundary_facet_axes() const { return m_boundary_facet_axes; }
    int boundary_facet_axis(size_t idx) const { return m_boundary_facet_axes.at(idx); }
    const coord_type &cell_coord(size_t idx) const { return m_cell_coords.at(idx); }
    const std::vector<coord_type> &cell_coords() const { return m_cell_coords; }

    int get_face_axis(int face_index) const;

  private:
    GridDatai m_cell_indices;// acts as a hash map for cell indexing based off a grid
    std::vector<coord_type> m_cell_coords;// reporst the cell found in each index (inverts m_cell_indices)
    std::vector<int> m_boundary_facet_axes;// parallel with DomainBoundary::boundary_facet_pairs
};


}// namespace mandoline
#include "mandoline/exterior_grid_impl.hpp"
