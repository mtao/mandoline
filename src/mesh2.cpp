#include "mandoline/mesh2.hpp"
#include <mtao/iterator/enumerate.hpp>
#include <mtao/logging/logger.hpp>
#include <mtao/geometry/mesh/halfedge.hpp>
#include <mtao/geometry/prune_vertices.hpp>
#include <mtao/geometry/kdtree.hpp>
#include <mtao/iterator/range.hpp>
#include "mandoline/diffgeo_utils.hpp"
#include "mandoline/line.hpp"
#include "mandoline/operators/boundary2.hpp"
#include "mandoline/operators/volume2.hpp"

namespace mandoline {
int CutCellMesh<2>::num_cutcells() const {
    return num_cutfaces();
}
int CutCellMesh<2>::num_cells() const {
    return num_faces();
}

int CutCellMesh<2>::num_cutfaces() const {
    return m_faces.size();
}
int CutCellMesh<2>::num_faces() const {
    return exterior_grid.num_cells() + num_cutfaces();
}
int CutCellMesh<2>::num_cutedges() const {
    return cut_edges().size();
}
int CutCellMesh<2>::num_edges() const {
    return exterior_grid.num_boundary_facets() + num_cutedges();
}

auto CutCellMesh<2>::centroids() const -> ColVecs {
    ColVecs C(2,num_cells());

    auto V = vertices();
    for(auto&& [idx, face]: mtao::iterator::enumerate(cut_faces())) {
        C.col(idx) = face.brep_centroid(V);
    }

    auto EC = C.rightCols(exterior_grid.num_cells());
    for(auto&& [idx,c]: mtao::iterator::enumerate(exterior_grid.cell_coords())) {
        EC.col(idx) = exterior_grid.cell_grid().vertex(c);
    }
    return C;

}

mtao::ColVectors<int, 3> CutCellMesh<2>::faces() const {
    std::vector<std::vector<int>> mycells;
    for (int i = 0; i < num_cells(); ++i) {
        for (auto &&c : cell(i)) {
            mycells.push_back(c);
        }
    }
    return mtao::geometry::mesh::earclipping(vertices(), mycells);
}



mtao::VecXi CutCellMesh<2>::regions() const {
    mtao::VecXi R(num_cells());
    for(auto&& [cid,c]: mtao::iterator::enumerate(cut_faces())) {
        R(cid) = c.region;
    }

    R.tail(exterior_grid.num_cells()) = mtao::eigen::stl2eigen(exterior_grid.regions());
    return R;
    
}

auto CutCellMesh<2>::volumes() const -> VecX {
    return operators::face_volumes(*this);
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


//auto CutCellMesh<2>::dual_vertices() const -> ColVecs {
//    ColVecs V = ColVecs::Zero(D, num_cells());
//    V.leftCols(StaggeredGrid::cell_size()) = StaggeredGrid::cell_vertices();
//
//    for (auto &&[i, c] : mtao::iterator::enumerate(hem.cell_halfedges())) {
//        int cell_index = hem.cell_index(c);
//        if (cell_index < StaggeredGrid::cell_size()) continue;
//
//        auto v = V.col(cell_index);
//        int count = 0;
//        mtao::geometry::mesh::cell_iterator(&hem, c)([&](auto &&e) {
//            auto u = vertex(e.vertex());
//            v += u;
//            ++count;
//        });
//        if (count > 0) {
//            v /= count;
//        }
//    }
//    return V;
//}
bool CutCellMesh<2>::in_cell(const VecCRef &p, int idx) const {
    auto V = vertices();
    return in_cell(V, p, idx);
}
bool CutCellMesh<2>::in_cell(const ColVecs &V, const VecCRef &p, int idx) const {
    if(idx < m_faces.size()) {
        return m_faces.at(idx).is_inside(V,p);
    } else {
        int c = exterior_grid.cell_index(p);
        if(c >= 0) {
            return true;
        }
    }
    return false;
}
int CutCellMesh<2>::cell_index(const VecCRef &p) const {
    auto [c, q] = StaggeredGrid::coord(p);
    if (!StaggeredGrid::cell_grid().valid_index(c)) {
        return -2;
    }


    auto V = vertices();
    int exterior_cell_index = exterior_grid.cell_index(c);
    if(exterior_cell_index == -1) {
        int grid_cell = grid_cell_index(c);

        if (auto it = cell_grid_ownership.find(grid_cell); it != cell_grid_ownership.end()) {
            auto &[i, cs] = *it;
            for (auto &&c : cs) {
                if (in_cell(V, p, c)) {
                    return c;
                }
            }
        }
    } else {
        return exterior_cell_index + num_cutfaces();
    }
    return -1;
}

std::set<std::vector<int>> CutCellMesh<2>::cell(int index) const {
    if (index < num_cutfaces()) {
        return m_faces.at(index).indices;
    } else if (int nidx = index - num_cutfaces(); nidx < exterior_grid.num_cells()) {
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
    /*

    int retind = -1;
    auto distance = [&](auto &&e) {

        auto vi = vertex(e.indices[0]);
        auto vj = vertex(e.indices[1]);
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

    if (!is_cutface_index(cind)) {

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
    */
    return {};
}

Eigen::SparseMatrix<double> CutCellMesh<2>::boundary(bool include_domain_boundary_faces) const {
    return operators::boundary(*this,include_domain_boundary_faces);


}
}// namespace mandoline
