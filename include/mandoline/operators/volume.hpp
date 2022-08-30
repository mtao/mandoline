#pragma once
#include "mandoline/exterior_grid.hpp"


namespace mandoline::operators {
template<int D>
balsa::eigen::VecXd dual_edge_lengths(const ExteriorGrid<D> &eg);
template<int D>
balsa::eigen::VecXd boundary_facet_volumes(const ExteriorGrid<D> &eg);
template<int D>
balsa::eigen::VecXd cell_volumes(const ExteriorGrid<D> &eg);

}// namespace mandoline::operators
