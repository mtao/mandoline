#pragma once
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mtao/geometry/grid/grid_data.hpp>
#include "mandoline/adaptive_grid.hpp"
#include <set>
#include "cutmesh.pb.h"
#include <balsa/eigen/stl2eigen.hpp>
#include <optional>
namespace mandoline::construction {


class AdaptiveGridFactory {
  public:
    constexpr static int logwidth = 1;
    constexpr static int width = 1 << logwidth;
    using Edge = std::array<int, 2>;
    using GridData3i = mtao::geometry::grid::GridDataD<int, 3>;
    using ActiveMask = std::bitset<1 << logwidth * 3>;
    using GridData3 = mtao::geometry::grid::GridDataD<bool, 3>;
    using GridData3b = mtao::geometry::grid::GridDataD<ActiveMask, 3>;
    using coord_type = std::array<int, 3>;
    using Cell = AdaptiveGrid::Cell;
    using AxialBEdgeMap = std::array<std::map<int, std::set<Edge>>, 3>;

    //cell mask
    AdaptiveGridFactory(const GridData3 &mask);


    std::tuple<std::array<std::set<Edge>, 3>, AxialBEdgeMap> compute_axial_edges(const std::optional<int> &max_level = {}) const;
    std::tuple<std::set<Edge>, AxialBEdgeMap> compute_edges(const std::optional<int> &max_level = {}) const;

    void make_cells(const std::optional<int> &max_level = {});

    using Indexer = mtao::geometry::grid::indexing::OrderedIndexer<3>;

    const static Indexer cmask_indexer;
    const static Indexer vmask_indexer;
    const static std::array<Indexer, 3> mask_edge_indexers;
    const static std::array<coord_type, 3> mask_edge_shapes;

    Indexer vertex_indexer;

    static int cmask_index(const coord_type &c) { return cmask_indexer.index(c); }
    static int cmask_index(int a, int b, int c) { return cmask_indexer.index(a, b, c); }
    static int vmask_index(const coord_type &c) { return vmask_indexer.index(c); }
    static int vmask_index(int a, int b, int c) { return vmask_indexer.index(a, b, c); }
    static int emask_index(int dim, const coord_type &c) { return mask_edge_indexers[dim].index(c); }
    static int emask_index(int dim, int a, int b, int c) { return mask_edge_indexers[dim].index(a, b, c); }
    int vertex_index(const coord_type &c) const { return vertex_indexer.index(c); }
    int vertex_index(int a, int b, int c) const { return vertex_indexer.index(a, b, c); }


    void make_edges(const std::optional<int> &max_level = {});
    void make_cells(const ActiveMask &mask, int level, const coord_type &coord);
    void make_cells(const GridData3 &mask, int level);

    void add_cell(const coord_type &corner, int jump);

    GridData3i grid_from_cells(const std::map<int, Cell> &cells) const;

    std::tuple<std::array<std::set<Edge>, 3>, AxialBEdgeMap> get_edges(const std::array<GridData3, 3> &edge_masks, int level, const coord_type &offset = {}) const;
    Edge get_edge(const coord_type &start, int jump, int dim) const;

    int get_jump(int level) const { return 1 << (level * logwidth); }
    int get_parent_jump(int level) const { return 1 << ((level + 1) * logwidth); }

    std::array<coord_type, 3> make_edge_shapes(const coord_type &coord) const;

    std::array<GridData3, 3> empty_edge_masks(const coord_type &shape, bool value = false) const;
    std::array<GridData3, 3> empty_edge_masks(int level, bool value = false) const;

    std::array<GridData3, 3> make_edge_mask(const GridData3 &mask, int level) const;
    std::array<GridData3, 3> make_edge_mask(const GridData3b &mask, int level) const;
    std::tuple<std::array<std::set<Edge>, 3>, AxialBEdgeMap> make_edges(const GridData3 &mask, int level) const;
    std::tuple<std::array<std::set<Edge>, 3>, AxialBEdgeMap> make_edges(const GridData3b &mask, int level) const;


    AdaptiveGrid create() const;


  public:
    balsa::eigen::ColVecs2i edges;
    balsa::eigen::ColVecs2i boundary_edges;
    std::vector<GridData3b> levels;
    std::vector<GridData3> levels_mask;
    GridData3 original;
    std::map<int, Cell> cells;

  private:
    template<typename GridAccessor>
    void set_edge_masks(const coord_type &shape, const coord_type &corner, std::array<GridData3, 3> &edge_mask, const GridAccessor &accessor) const;
};

}// namespace mandoline::construction
