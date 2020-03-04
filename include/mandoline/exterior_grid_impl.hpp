#pragma once
#include "mandoline/exterior_grid.hpp"


namespace mandoline {
template<int D>
ExteriorGrid<D>::ExteriorGrid(const Base &g, const GridDatab &cell_mask) : Base(g), m_cell_indices(cell_mask.shape()) {
    assert(g.cell_shape() == cell_mask.shape());
    using namespace mtao::geometry::grid;
    int counter = 0;
    std::transform(cell_mask.begin(), cell_mask.end(), m_cell_indices.begin(), [&counter](bool outside) -> int {
        if (!outside) {
            return -1;
        } else {
            return counter++;
        }
    });

    m_cell_coords.resize(counter);
    m_cell_indices.loop([&](auto &&c, int idx) {
        if (idx >= 0) {
            m_cell_coords[idx] = c;
        }
    });

    // we'll assume that the grid is mostly full
    m_boundary_facet_pairs.reserve((D - 1) * cell_mask.size());
    m_boundary_facet_axes.reserve((D - 1) * cell_mask.size());
    for (int i = 0; i < D; ++i) {
        coord_type s = cell_shape();
        s[i] = 1;
        // goes from -1 to 0
        utils::multi_loop(s, [&](auto v) {
            int pi = m_cell_indices(v);
            if (pi >= 0) {
                m_boundary_facet_pairs.emplace_back(std::array<int, 2>{ { -2, pi } });
                m_boundary_facet_axes.emplace_back(i);
            }
        });
        s = cell_shape();
        s[i] -= 1;
        if (s[i] < 1) continue;
        // goes from 0 to s[i]-1
        utils::multi_loop(s, [&](auto v) {
            auto v2 = v;
            v2[i] += 1;
            int ni = m_cell_indices(v);
            int pi = m_cell_indices(v2);
            if (ni == -1 || pi == -1) {
                return;
            } else {
                m_boundary_facet_pairs.emplace_back(std::array<int, 2>{ { ni, pi } });
                m_boundary_facet_axes.emplace_back(i);
            }
        });
        s[i] = 1;
        // goes from s[i]-1 to s[i]
        utils::multi_loop(s, [&](auto v) {
            v[i] = cell_shape()[i] - 1;
            int ni = m_cell_indices(v);
            if (ni >= 0) {
                m_boundary_facet_pairs.emplace_back(std::array<int, 2>{ { ni, -2 } });
                //if constexpr (D == 2) {
                //    m_boundary_facet_axes.emplace_back(1 - i);
                //} else {
                m_boundary_facet_axes.emplace_back(i);
                //}
            }
        });
    }
}

template<int D>
std::vector<Eigen::Triplet<double>> ExteriorGrid<D>::boundary_facet_to_staggered_grid(int offset) const {
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(boundary_facet_pairs().size());
    for (auto &&[row, pr, axis] : mtao::iterator::enumerate(boundary_facet_pairs(), boundary_facet_axes())) {
        const int gaxis = (D == 2) ? 1 - axis : axis;
        auto [ai, bi] = pr;
        if (ai == -2) {
            auto b = cell_coord(bi);
            assert(b[axis] == 0);
            int col = Base::template staggered_index<D - 1>(b, gaxis);
            auto g = Base::template grid<1>(axis);
            trips.emplace_back(row + offset, col, 1);
        } else if (bi == -2) {
            auto a = cell_coord(ai);
            a[axis] += 1;
            assert(a[axis] <= vertex_shape()[axis]);
            int col = Base::template staggered_index<D - 1>(a, gaxis);
            auto g = Base::template grid<1>(axis);
            trips.emplace_back(row + offset, col, 1);
        } else if (ai >= 0 && bi >= 0) {
            auto a = cell_coord(ai);
            auto b = cell_coord(bi);
            assert(b[axis] - 1 == a[axis]);
            auto g = Base::template grid<1>(axis);
            int col = Base::template staggered_index<D - 1>(b, gaxis);
            trips.emplace_back(row + offset, col, 1);
        }
    }
    return trips;
}
template<int D>
mtao::VecXd ExteriorGrid<D>::boundary_facet_volumes(bool make_boundary) const {
    auto &A = boundary_facet_axes();
    mtao::VecXd ret(A.size());
    ret.setZero();
    auto vols = Base::template form_volumes<D - 0>();

    for (auto &&[row, axis] : mtao::iterator::enumerate(boundary_facet_axes())) {
        if (make_boundary || !is_boundary_facet(row)) {
            ret(row) = vols[axis];
        }
    }
    return ret;
}

}// namespace mandoline
