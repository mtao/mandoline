#pragma once
#include "mandoline/operators/volume.hpp"


namespace mandoline::operators {
template<int D>
balsa::eigen::VecXd dual_edge_lengths(const ExteriorGrid<D> &eg) {
    balsa::eigen::VecXd R(eg.boundary_facet_size());
    using Base = typename ExteriorGrid<D>::Base;
    auto grid_vols = eg.Base::template form_volumes<1>();
    for (auto &&[v, axis] : mtao::iterator::zip(mtao::iterator::shell(R.data(), R.data() + R.size()), eg.boundary_facet_pairs())) {
        v = grid_vols[axis];
    }
    return R;
}

template<int D>
balsa::eigen::VecXd boundary_facet_volumes(const ExteriorGrid<D> &eg) {
    balsa::eigen::VecXd R(eg.boundary_facet_size());
    using Base = typename ExteriorGrid<D>::Base;
    auto grid_vols = eg.Base::template form_volumes<D - 1>();
    // in 2d we have this issue wehre
    if constexpr (D == 2) {
        for (auto &&[v, axis] : mtao::iterator::zip(mtao::iterator::shell(R.data(), R.data() + R.size()), eg.boundary_facet_axes())) {
            v = grid_vols[1 - axis];
        }
    } else {
        for (auto &&[v, axis] : mtao::iterator::zip(mtao::iterator::shell(R.data(), R.data() + R.size()), eg.boundary_facet_pairs())) {
            v = grid_vols(axis);
        }
    }
    return R;
}
template<int D>
balsa::eigen::VecXd cell_volumes(const ExteriorGrid<D> &eg) {
    auto vol = eg.dx().prod();
    return balsa::eigen::VecXd::Constant(eg.num_cells(),vol);
}
}// namespace mandoline::operators
