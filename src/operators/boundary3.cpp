#include "mandoline/operators/boundary3.hpp"
#include "mandoline/operators/boundary.hpp"


namespace mandoline::operators {

Eigen::SparseMatrix<double> boundary(const CutCellMesh<3> &ccm, bool include_domain_boundary) {
    auto trips = boundary_triplets(ccm.adaptive_grid(), ccm.faces().size(), include_domain_boundary);
    Eigen::SparseMatrix<double> B(ccm.face_size(), ccm.cell_size());

    auto g = ccm.adaptive_grid().cell_ownership_grid();

    for (auto &&c : ccm.cells()) {
        int region = c.region;
        for (auto &&[fidx, s] : c) {
            auto &f = ccm.faces()[fidx];
            if(include_domain_boundary || !(f.external_boundary && std::get<0>(*f.external_boundary) == -2)) {
                trips.emplace_back(fidx, c.index, s ? 1 : -1);
            }
        }
    }
    for (auto &&[fidx, f] : mtao::iterator::enumerate(ccm.faces())) {
        if (f.external_boundary) {
            auto [c, s] = *f.external_boundary;
            if (c >= 0) {
                trips.emplace_back(fidx, g.get(c), s ? -1 : 1);
            }
        }
    }


    B.setFromTriplets(trips.begin(), trips.end());
    return B;
}
//for removing mesh faces
std::set<int> grid_boundary_faces(const CutCellMesh<3> &ccm) {
    std::set<int> ret;
    for (auto &&[fidx, f] : mtao::iterator::enumerate(ccm.faces())) {
        auto mask = f.mask();
        for (auto [dim, valo] : mtao::iterator::enumerate(mask)) {
            if (valo) {
                int val = *valo;
                if (val == 0) {
                    ret.insert(fidx);
                } else if (val == ccm.vertex_shape()[dim] - 1) {
                    ret.insert(fidx);
                }
            }
        }
    }
    // adaptive grid part
    int fidx_offset = ccm.cut_face_size();
    auto adret = grid_boundary_faces(ccm.adaptive_grid(),fidx_offset);
    ret.merge(adret);
    return ret;
}
std::set<int> grid_boundary_faces(const AdaptiveGrid& db, int offset) {
    std::set<int> ret;
    for (auto &&[fidx_, f] : mtao::iterator::enumerate(db.faces())) {
        int fidx = fidx_ + offset;
        int dim = f.dimension();
        int val = f.corner()[dim];
        if (val == 0) {
            ret.insert(fidx);
        } else if (val == db.cell_shape()[dim]) {
            ret.insert(fidx);
        }
    }
    return ret;
}

std::vector<Eigen::Triplet<double>> boundary_triplets(const AdaptiveGrid& ag, int offset, bool domain_boundary ) {
    std::vector<Eigen::Triplet<double>> trips;
    for (auto &&[i, face] : mtao::iterator::enumerate(ag.faces())) {
        auto &e = face.dual_edge;
        if(!domain_boundary || ag.is_boundary_face(e)) {
            auto [l, h] = face.dual_edge;
            if (l >= 0) {
                trips.emplace_back(offset + i, l, -1);
            }
            if (h >= 0) {
                trips.emplace_back(offset + i, h, 1);
            }
        }
    }
    return trips;
}
}// namespace mandoline::operators
