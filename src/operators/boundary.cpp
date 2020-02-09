#include "mandoline/operators/boundary.hpp"


namespace mandoline::operators {

    Eigen::SparseMatrix<double> boundary(const CutCellMesh<3>& ccm) {
        auto trips = ccm.adaptive_grid().boundary_triplets(ccm.faces().size());
        Eigen::SparseMatrix<double> B(ccm.face_size(),ccm.cell_size());

        auto g = ccm.adaptive_grid().cell_ownership_grid();

        for(auto&& c: ccm.cells()) {
            int region = c.region;
            for(auto&& [fidx,s]: c) {
                auto& f = ccm.faces()[fidx];
                if(f.is_axial_face()) {
                    int a = ccm.cell_shape()[f.as_axial_axis()];
                    int v = f.as_axial_coord();
                    if(v > 0 && v < a) {
                        trips.emplace_back(fidx,c.index, s?1:-1);
                    }
                } else {
                    trips.emplace_back(fidx,c.index, s?1:-1);
                }
            }
        }
        for(auto&& [fidx,f]: mtao::iterator::enumerate(ccm.faces())) {
            if(f.external_boundary) {
                auto [c,s] = *f.external_boundary;
                if(c >= 0) {
                    trips.emplace_back(fidx,g.get(c), s?-1:1);
                }
            }
        }


        B.setFromTriplets(trips.begin(),trips.end());
        return B;
    }
    //for removing mesh faces
    std::set<int> grid_boundary_faces(const CutCellMesh<3>& ccm) {
        std::set<int> ret;
        std::cout << "Going through faces" << std::endl;
        for(auto&& [fidx,f]: mtao::iterator::enumerate(ccm.faces())) {
            auto mask = f.mask();
            for(auto [dim,valo]: mtao::iterator::enumerate(mask)) {
                if(valo) {
                    int val = *valo;
                    if(val == 0) {
                        ret.insert(fidx);
                    } else if(val == ccm.vertex_shape()[dim]-1) {
                        ret.insert(fidx);
                    }
                }
            }
        }
        // adaptive grid part
        int fidx_offset = ccm.cut_face_size();
        auto adret = ccm.adaptive_grid().grid_boundary_faces(fidx_offset);
        ret.merge(adret);
        return ret;
    }
}
