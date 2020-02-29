#include "mandoline/mesh2.hpp"
#include "mandoline/operators/boundary.hpp"


namespace mandoline::operators {

Eigen::SparseMatrix<double> boundary(const CutCellMesh<2> &ccm, bool include_domain_boundary) {
    auto triplets = ccm.exterior_grid.boundary_triplets(include_domain_boundary);
    int num_cedges = ccm.cut_edges().size();
    int num_cfaces = ccm.cut_faces().size();
    std::transform(triplets.begin(),triplets.end(), triplets.begin(), [&](const Eigen::Triplet<double>& trip) {
            return Eigen::Triplet<double>{num_cedges+trip.row(), num_cfaces+trip.col(),trip.value()};
            });

    for(auto&& [fidx, emap]: ccm.m_face_boundary_map) {

        for(auto&& [eidx,sgn]: emap) {
            auto& e = ccm.cut_edges()[eidx];
            // if we have include_domain condition then we DO NOT want boundaries taht connect to nothing
            if(include_domain_boundary || !(e.external_boundary && std::get<0>(*e.external_boundary) == -2)) {
                triplets.emplace_back(eidx,fidx,sgn?-1:1);
            }
        }
    }
    for(auto&& [idx,e]: mtao::iterator::enumerate(ccm.cut_edges())) {
        if(e.external_boundary) {
            auto [of,sgn] = *e.external_boundary;
            if(of >= 0) {
                triplets.emplace_back(idx,ccm.exterior_grid.cell_indices().get(of)+num_cfaces,sgn?-1:1);
            }
        }
    }

    Eigen::SparseMatrix<double>  A(ccm.num_edges(),ccm.num_faces());
    A.setFromTriplets(triplets.begin(),triplets.end());
    return A;
}

std::set<int> grid_boundary_faces(const CutCellMesh<2> &ccm) {
    std::set<int> ret;
    for (auto &&[fidx, f] : mtao::iterator::enumerate(ccm.cut_edges())) {
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
    int fidx_offset = ccm.cut_edge_size();
    auto adret = grid_boundary_faces(ccm.exterior_grid, fidx_offset);

    ret.merge(adret);
    return ret;
}

}
