#include "mandoline/tools/exploded_mesh.hpp"

namespace mandoline::tools {
    MeshExploder::MeshExploder(const CutCellMesh<3>& ccm) {
        auto V = ccm.vertices();
        int off = 0;
        for(auto&& [i,c]: mtao::iterator::enumerate(ccm.cells())) {
            auto [v,f] = c.get_mesh(V,ccm.faces());
            Vs.emplace_back(std::move(v));
            Fs.emplace_back(std::move(f));
            Cs.emplace_back(c.centroid(V,ccm.faces()));
        }
    }

    std::vector<int> MeshExploder::offsets() const {
        std::vector<int> ret;
        ret.reserve(Vs.size()+1);
        int off = 0;
        ret.push_back(off);
        for(auto&& v: Vs) {
            off += v.cols();
            ret.push_back(off);
        }
        return ret;
    }

    std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> MeshExploder::mesh(double scale) const {
       return {V(scale),F()};
    }
    mtao::ColVecs3d MeshExploder::V(double scale) const {

        std::vector<mtao::ColVecs3d> V2(Vs.size());
        for(auto&& [V,C]: mtao::iterator::zip(Vs,Cs)) {
            V2.push_back(V.colwise() + (scale-1) * C);
        }
        return mtao::eigen::hstack_iter(V2.begin(),V2.end());
    }
    mtao::ColVecs3i MeshExploder::F() const {

        std::vector<mtao::ColVecs3i> F2(Fs.size());
        for(auto&& [F,o]: mtao::iterator::zip(Fs,offsets())) {
            F2.push_back(F.array() + o);
        }
        return mtao::eigen::hstack_iter(F2.begin(),F2.end());
    }
}
