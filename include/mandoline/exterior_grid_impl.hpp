#pragma once
#include "mandoline/exterior_grid.hpp"


namespace mandoline {
    template <int D>
        ExteriorGrid<D>::ExteriorGrid(const GridDatab& cell_mask): m_cell_indices(cell_mask.shape()) {
            using namespace mtao::geometry::grid;
            int counter = 0;
            std::transform(cell_mask.begin(),cell_mask.end(), m_cell_indices.begin(),[&counter](bool inside) -> int {
                    if(inside) {
                    return -1;
                    } else {
                    return counter++;
                    }
                    });

            m_cell_coords.reserve(counter);
            m_cell_indices.loop([&](auto&& c,auto) {
                    m_cell_coords.emplace_back(c);
                    });

            // we'll assume that the grid is mostly full
            m_boundary_facet_pairs.reserve((D-1) * cell_mask.size());
            m_boundary_facet_axes.reserve((D-1) * cell_mask.size());
            for(int i = 0; i < D; ++i) {
                coord_type s = cell_shape();
                s[i]=1;
                utils::multi_loop(s,[&](auto v) {
                        int pi = m_cell_indices(v);
                        if(pi >= 0) {
                        m_boundary_facet_pairs.emplace_back(std::array<int,2>{{-2,pi}});
                        m_boundary_facet_axes.emplace_back(i);
                        }
                        });
                s = cell_shape();
                s[i]-=2;
                if(s[i] < 1) continue;
                utils::multi_loop(s,[&](auto v) {
                        auto v2 = v;
                        v[i]+=1;
                        v2[i]+=2;
                        int ni = m_cell_indices(v);
                        int pi = m_cell_indices(v2);
                        if(ni == -1 && pi == -1) {
                        return;
                        } else {
                        m_boundary_facet_pairs.emplace_back(std::array<int,2>{{ni,pi}});
                        m_boundary_facet_axes.emplace_back(i);
                        }
                        });
                s[i]=1;
                utils::multi_loop(s,[&](auto v) {
                        v[i] = cell_shape()[i]-1;
                        int ni = m_cell_indices(v);
                        if(ni >= 0) {
                        m_boundary_facet_pairs.emplace_back(std::array<int,2>{{ni,-2}});
                        m_boundary_facet_axes.emplace_back(i);
                        }
                        });
            }

        }
    template <int D>
        mtao::VecXd ExteriorGrid<D>::face_volumes(bool mask_boundary) const {
            mtao::VecXd R(num_faces());
            auto&dx = Base::dx();
            mtao::Vec3d dws = mtao::Vec3d::Ones();
            for(int i = 0; i < D; ++i) {
                for(int j = 0; j < D-1; ++j) {
                    dws(i) *= dx((j+i)%D);
                }
            }


            R.setZero();
            for(int i = 0; i < R.size(); ++i) {
                if(!(mask_boundary && is_boundary_facet(i))) {
                    R(i) = dws(boundary_facet_axis(i));
                }
            }
            return R;
        }


}
