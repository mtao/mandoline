#include "mandoline/mesh2.hpp"
#include <mtao/iterator/enumerate.hpp>
#include <mtao/logging/logger.hpp>
#include <mtao/geometry/mesh/halfedge.hpp>
#include <mtao/geometry/prune_vertices.hpp>
#include <mtao/geometry/kdtree.hpp>
#include <mtao/iterator/range.hpp>
#include "mandoline/diffgeo_utils.hpp"
#include "mandoline/line.hpp"

namespace mandoline {
int CutCellMesh<2>::num_cutcells() const {
    return m_faces.size();
}
int CutCellMesh<2>::num_cells() const {
    return exterior_grid.num_cells() + num_cutcells();
}

mtao::ColVectors<int, 3> CutCellMesh<2>::faces() const {
    std::cout << "True cutfaces: " << std::endl;
    for (auto &&[i, c] : mtao::iterator::enumerate(m_faces)) {
        std::cout << i << ")): ";
        for (auto &&v : c.indices) {
            for (auto &&v : v) {
                std::cout << v << ",";
            }
            std::cout << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "HEM cutfaces: " << std::endl;
    for (auto &&[i, c] : mtao::iterator::enumerate(hem.cells())) {
        std::cout << i << ")): ";
        for (auto &&v : c) {
            std::cout << v << ",";
        }
        std::cout << std::endl;
    }
    std::vector<std::vector<int>> mycells;
    for (int i = 0; i < num_cells(); ++i) {
        for (auto &&c : cell(i)) {
            mycells.push_back(c);
        }
    }
    return mtao::geometry::mesh::earclipping(vertices(), mycells);
}
auto CutCellMesh<2>::volumes() const -> VecX {
    VecX V(num_cells());
    V.topRows(StaggeredGrid::cell_size()).array() = dx().prod();

    auto C = hem.cell_halfedges();

    for (auto &&c : hem.cell_halfedges()) {
        int cell_index = hem.cell_index(c);
        if (is_grid_face(cell_index)) continue;

        double &a = V(cell_index) = 0;
        mtao::geometry::mesh::cell_iterator(&hem, c)([&](auto &&e) {
            auto u = vertex(e.vertex());
            auto v = vertex(e.get_dual().vertex());
            a += .5 * (u[0] * v[1] - v[0] * u[1]);
        });
        a = std::abs(a);
    }
    return V;
}
auto CutCellMesh<2>::dual_edge(int idx) const -> Edge {
    if (is_grid_edge(idx)) {
        return grid_dual_edge(idx);
    } else {
        auto he = halfedges_per_edge.col(idx);
        return { { hem.cell_index(he(0)), hem.cell_index(he(1)) } };
    }
}
auto CutCellMesh<2>::dual_edge_volumes() const -> VecX {
    VecX V(edge_size());

    auto offsets = StaggeredGrid::template offsets<1>();
    auto mdx = dx();
    for (int i = 0; i < D; ++i) {
        double odx = 1;
        for (int j = 0; j < D - 1; ++j) {
            odx *= mdx((i + j + 1) % D);
        }
        V.segment(offsets[i], StaggeredGrid::template staggered_size<1>(i)).array() = odx;
    }

    return V;
}


auto CutCellMesh<2>::dual_vertices() const -> ColVecs {
    ColVecs V = ColVecs::Zero(D, num_cells());
    V.leftCols(StaggeredGrid::cell_size()) = StaggeredGrid::cell_vertices();

    for (auto &&[i, c] : mtao::iterator::enumerate(hem.cell_halfedges())) {
        int cell_index = hem.cell_index(c);
        if (cell_index < StaggeredGrid::cell_size()) continue;

        auto v = V.col(cell_index);
        int count = 0;
        mtao::geometry::mesh::cell_iterator(&hem, c)([&](auto &&e) {
            auto u = vertex(e.vertex());
            v += u;
            ++count;
        });
        if (count > 0) {
            v /= count;
        }
    }
    return V;
}
bool CutCellMesh<2>::in_cell(const VecCRef &p, int idx) const {
    auto V = vertices();
    return in_cell(V, p, idx);
}
bool CutCellMesh<2>::in_cell(const ColVecs &V, const VecCRef &p, int idx) const {
    double wn = 0;
    for (auto &&f : m_faces.at(idx).indices) {
        wn += mtao::geometry::interior_winding_number(V, f, p);
    }
    return std::abs(wn) > .5;
}
int CutCellMesh<2>::cell_index(const VecCRef &p) const {
    auto [c, q] = StaggeredGrid::coord(p);
    if (!StaggeredGrid::cell_grid().valid_index(c)) {
        return -1;
    }

    auto V = vertices();
    int grid_cell = StaggeredGrid::cell_index(c);

    if (auto it = cell_grid_ownership.find(grid_cell); it != cell_grid_ownership.end()) {
        auto &[i, cs] = *it;
        for (auto &&c : cs) {
            if (in_cell(V, p, c)) {
                int ret = c;
                return ret;
            }
        }
    }
    return grid_cell;
}

std::set<std::vector<int>> CutCellMesh<2>::cell(int index) const {
    if (index < num_cutcells()) {
        return m_faces.at(index).indices;
    } else if (int nidx = index - num_cutcells(); nidx < exterior_grid.num_cells()) {
        std::vector<int> r;

        coord_type cidx = exterior_grid.cell_coord(nidx);
        r.push_back(vertex_index(std::array<int, 2>{ { cidx[0] + 0, cidx[1] + 0 } }));
        r.push_back(vertex_index(std::array<int, 2>{ { cidx[0] + 0, cidx[1] + 1 } }));
        r.push_back(vertex_index(std::array<int, 2>{ { cidx[0] + 1, cidx[1] + 1 } }));
        r.push_back(vertex_index(std::array<int, 2>{ { cidx[0] + 1, cidx[1] + 0 } }));
        return { r };

    } else {
        std::cout << index << "/" << num_cutcells() << "/" << exterior_grid.num_cells() << std::endl;
        return {};
    }
}

int CutCellMesh<2>::nearest_edge_index(const VecCRef &p) const {
    int cind = cell_index(p);

    int retind = -1;
    auto distance = [&](auto &&e) {
        int i = e.vertex();
        int j = e.get_dual().vertex();
        auto vi = vertex(i);
        auto vj = vertex(j);
        Line<2> l(vi, vj);
        return (l.project(p) - p).cwiseQuotient(dx()).norm();
    };
    double mindist = std::numeric_limits<double>::max();
    auto update_min = [&](double d, int ind) {
        if (d < mindist) {
            mindist = d;
            retind = ind;
        }
    };

    if (is_grid_face(cind)) {

        auto [c, q] = coord(p);
        auto &&[i, j] = c;
        auto &&[a, b] = q;

        auto active_predicate = [&](int i, int j) {
            return !active_grid_cell_mask().valid_index(i, j) || active_grid_cell_mask()(i, j);
        };
        if (active_predicate(i, j)) {
            //b / y
            if (active_predicate(i, j + 1)) {
                update_min(1 - b, v_index({ { i, j + 1 } }));
            }
            if (active_predicate(i, j - 1)) {
                update_min(b, v_index(c));
            }

            //a /xy
            if (active_predicate(i + 1, j)) {
                update_min(1 - a, u_index({ { i + 1, j } }));
            }
            if (active_predicate(i - 1, j)) {
                update_min(a, u_index(c));
            }
        }
    }
    auto es = hem.cell_edges(cind);
    if (!es.empty()) {

        for (auto &&e : es) {
            update_min(distance(e), halfedge_to_edge_index.at(int(e)) + StaggeredGrid::edge_size());
        }
    }
    return retind;
}
}// namespace mandoline
