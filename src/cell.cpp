#include "mandoline/cell.hpp"
#include <mtao/logging/logger.hpp>
#include <mtao/eigen/stack.h>

namespace mandoline {
    std::vector<Eigen::Triplet<double>> CutCell::boundary_triplets() const {
        std::vector<Eigen::Triplet<double>> trips;
        for(auto&& [i,o]: *this) {
            trips.emplace_back(i,index,o?1:-1);
        }
        return trips;

    }
    /*
       auto CutCell::grid_cell(const std::vector<CutFace<3>>& F) const -> std::array<int,3> {
       if(F.empty()) { return {}; }
       if(empty()) { return {}; }
       auto mask = [&](int idx) {
       return F[idx].mask();
       };
       auto possibles = [](const coord_mask& c) -> std::set<coord_type> {

       std::set<coord_type> ret;
       for(int i = 0; i < (2 << 3); ++i) {
       std::bitset<D> bs(i);
       if((bs&clamped_indices) == bs) {

       coord_type cc = c;
       for(int j = 0; j < D; ++j) {
       cc[j] -= bs[j]?1:0;
       }
       ret.emplace(cc);

       }
       }
       return ret;
       };
       std::set<coord_type> possibles = possibles(std::get<0>(*this.begin()));

       for(auto&& f: F) {
       auto m = 
       auto s = GV(f).possible_cells();
       std::set<CoordType> i;
       std::set_intersection(possibles.begin(),possibles.end(),s.begin(),s.end(),std::inserter(i,i.end()));
       possibles = std::move(i);
       if(possibles.empty()) {
       return {};
       }
       }

       return possibles;
       }
       */
    void  CutCell::serialize(CutMeshProto::Cell& cell) const {
        cell.set_id(index);
        cell.set_region(region);
        protobuf::serialize(grid_cell,*cell.mutable_grid_cell());
        auto&& pmap = *cell.mutable_entries();
        for(auto&& [a,b]: *this) {
            pmap[a] = b;
        }
    }
    CutCell CutCell::from_proto(const CutMeshProto::Cell& cell) {
        CutCell c;
        auto&& entries = cell.entries();
        //c = std::map<int,bool>(entries.begin(),entries.end());
        for(auto&& [a,b]: entries) {
            c[a] = b;
        }
        c.index = cell.id();
        c.region = cell.region();
        protobuf::deserialize(cell.grid_cell(),c.grid_cell);
        return c;
    }
    mtao::ColVecs3i CutCell::triangulated(const std::vector<CutFace<3>>& Fs) const {
        std::vector<mtao::ColVecs3i> T;
        for(auto&& [s,b]: *this) {
            auto&& F = Fs[s];
            if(F.triangulation) {
                auto&& t = *F.triangulation;
                T.emplace_back(t);
                if(!b) {
                    auto&& b = T.back();
                    auto tmp = b.row(0).eval();
                    b.row(0) = b.row(1);
                    b.row(1) = tmp;
                }
            } else {
                mtao::logging::error() << "Cell tried to fetch triangulation from an untriangulated face";
            }
        }

        return mtao::eigen::hstack_iter(T.begin(),T.end());
    }

    double CutCell::volume(const mtao::VecXd& face_brep_vols) const {

        double vol = 0;

        for(auto&& [f,b]: *this) {
            double sign = b?1:-1;
            vol += sign * face_brep_vols(f);

        }
        return vol;
    }
    mtao::Vec3d CutCell::moment(const mtao::ColVecs3d& face_brep_cents) const {

        mtao::Vec3d c = mtao::Vec3d::Zero();

        for(auto&& [f,b]: *this) {
            double sign = b?1:-1;
            //c += sign * face_brep_cents.col(f);
            c += face_brep_cents.col(f);

        }
        c /= size();
        return c;
    }
}
