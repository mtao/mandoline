#include "cutmesh_validation.hpp"

using namespace mandoline;

bool paired_boundary(const CutCellMesh<3>& ccm) {
    mtao::VecXd FC = mtao::VecXd::Zero(ccm.faces().size());
    mtao::VecXd FB = mtao::VecXd::Zero(ccm.faces().size());
    std::map<int,mtao::VecXd> region_counts;
    std::map<int,mtao::VecXd> region_boundaries;
    auto regions = ccm.regions();
    for(auto&& r: regions) {
        if(region_counts.find(r) == region_counts.end()) {
            region_counts[r] = mtao::VecXd::Zero(ccm.faces().size());
            region_boundaries[r] = mtao::VecXd::Zero(ccm.faces().size());
        }
    }
    for(auto&& c: ccm.cells()) {
        int region = c.region;
        auto& counts = region_counts[region];
        auto& boundaries = region_boundaries[region];
        for(auto&& [f,s]: c) {
            counts[f]++;
            boundaries[f] += s?1:-1;
        }
    }
    for(auto&& [r,counts]: region_counts) {
        for(int i = 0; i < counts.size(); ++i) {
            if(counts[i] == 1) {
                if(!ccm.is_mesh_face(i) && ccm.faces()[i].indices.begin()->size() != 4) {
                    return false;
                }
            } else if(counts[i] > 2) {
                return false;
            }
        }
        FC += counts;
    }

    for(auto&& [r,boundaries]: region_boundaries) {
        for(int i = 0; i < boundaries.size(); ++i) {
            if(std::abs(i) == 1) {
                if(!ccm.is_mesh_face(i) && ccm.faces()[i].indices.begin()->size() != 4) {
                    return false;
                }
            }
        }
        FB += boundaries;
    }
    for(int i =0 ; i< FB.rows(); ++i ) {
        if(FB(i) < -1 || FB(i) > 1 && FC(i) < 0 || FC(i) > 2) {
            return false;
        } else {
            using CoordType = std::array<int,3>;
            if(FC(i) == 1) {
                if(ccm.folded_faces().find(i) == ccm.folded_faces().end()) {
                    if(!bool(ccm.faces()[i].external_boundary)) {
                        for(auto&& inds: ccm.faces()[i].indices) {

                            auto pc = possible_cells(ccm,inds);
                            if(pc.size() != 2) {
                                continue;
                            }
                            std::array<CoordType,2> pca;
                            std::copy(pc.begin(),pc.end(),pca.begin());
                            //find the boundary cells axis
                            int idx;
                            for(idx=0; idx < 3; ++idx) {
                                if(pca[0][idx] != pca[1][idx]) {
                                    break;
                                }
                            }
                            int val = pca[1][idx];


                            if(val != 0 && val != ccm.vertex_shape()[idx]-1) {
                                return false;
                            }
                        }
                    }
                }
            }
        }

    }
    return true;
}
