#pragma once
#include "mandoline/cutface.hpp"

namespace mandoline {

    template <int D>
        bool CutFace<D>::is_mesh_face() const {
            return std::holds_alternative<int>(id);
        }
    template <int D>
        bool  CutFace<D>::is_axial_face() const {
            return std::holds_alternative<std::array<int,2>>(id);
        }
    template <int D>
        int  CutFace<D>::as_face_id() const {
            return std::get<int>(id);
        }
    template <int D>
    const std::array<int,2>&  CutFace<D>::as_axial_id() const {
        return std::get<std::array<int,2>>(id);
    }
    template <int D>
    int  CutFace<D>::as_axial_axis() const {
        return std::get<std::array<int,2>>(id)[0];
    }
    template <int D>
    int  CutFace<D>::as_axial_coord() const {
        return std::get<std::array<int,2>>(id)[1];
    }
    template <int D>
        CutFace<D>::operator std::string() const {
            std::stringstream ss;
            ss << "{";

            std::visit([&](auto&& f) {
                    using T = std::decay_t<decltype(f)>;
                    if constexpr(std::is_same_v<T,int>) {
                    ss << "(" << f << ")";
                    } else {
                    auto [a,b] = f;
                    ss << "(" << a << ":" << b << ")";

                    }
                    },id);
            for(auto&& f: indices) {
                ss << "[";
                std::copy(f.begin(),f.end(),std::ostream_iterator<int>(ss,","));
                ss << "],";
            }
            ss << "}";


            return ss.str();

        }

    template <int D>
        void   CutFace<D>::serialize(CutMeshProto::CutFace& face) const {
            protobuf::serialize(N,*face.mutable_normal());
            if(is_mesh_face()) {
                face.set_face_id(std::get<int>(id));
            } else {
                auto& ap = *face.mutable_plane_id();
                auto&& [a,b] = std::get<std::array<int,2>>(id);
                ap.set_axis(a);
                ap.set_value(b);
            }
            for(auto&& curve: indices) {
                auto&& c = *face.add_curves();
                for(auto&& v: curve) {
                    c.add_indices(v);
                }
            }
            if(triangulation) {
                auto&& T = *triangulation;
                for(int i = 0; i < T.cols(); ++i) {
                    protobuf::serialize(T.col(i),*face.add_triangulation());
                }
            }
            if(external_boundary) {
                auto [b,s] = *external_boundary;
                auto&& fb = *face.mutable_face_boundary();
                fb.set_index(b);
                fb.set_sign(s);
            }
        }

    template <int D>
        CutFace<D>  CutFace<D>::from_proto(const CutMeshProto::CutFace& face) {
            CutFace<D> ret;
            ret.N = protobuf::deserialize(face.normal());
            if(face.id_case() == CutMeshProto::CutFace::IdCase::kFaceId) {
                ret.id = face.face_id();
            } else {
                auto&& pid = face.plane_id();
                ret.id = std::array<int,2>{{pid.axis(),pid.value()}};
            }
            for(auto&& c: face.curves()) {
                std::vector<int> curve(c.indices().size());
                std::copy(c.indices().begin(),c.indices().end(),curve.begin());
                ret.indices.insert(curve);
            }

            if(face.triangulation().size() > 0) {
                mtao::ColVecs3i T;
                T.resize(3,face.triangulation().size());
                for(int i = 0; i < T.cols(); ++i) {
                    T.col(i) = protobuf::deserialize(face.triangulation(i));
                }
                ret.triangulation = T;
            }
            if(face.has_face_boundary()) {
                auto&& fb = face.face_boundary();
                ret.external_boundary = {fb.index(),fb.sign()};
            }
            return ret;
        }

    template <int D>
        template <typename Derived>
        mtao::Vec3d CutFace<D>::brep_centroid(const Eigen::MatrixBase<Derived>& V, bool use_triangulation) const {
            mtao::Vec3d ret = mtao::Vec3d::Zero();
            double tot_vol = 0;
            int count = 0;
            for(auto&& i: indices) {
                count += i.size();
                for(auto&& i: i) {
                    ret += V.col(i);
                }
                /*
                   auto vol = brep_volume(V,use_triangulation);
                   tot_vol += vol;
                   mtao::Vec3d cent = mtao::geometry::curve_centroid(V,i);
                   ret += vol * cent;
                   */
            }
            /*
               if(tot_vol > 1e-10) {
               ret /= tot_vol;
               }
               */
            ret /= count;
            return ret;
        }










    template <int D>
        template <typename Derived>
        double CutFace<D>::brep_volume(const Eigen::MatrixBase<Derived>& V, bool use_triangulation) const {

            if(use_triangulation) {
                return mtao::geometry::brep_volume(V,*triangulation);
            }
            double vol = 0;
            auto loop_volume = [&](const std::vector<int>& loop) {

                double vol = 0;
                mtao::Matrix<double,3,4> S;
                S.col(3).setZero();
                S.col(0) = V.col(loop[0]);
                S.col(2) = V.col(loop[1]);

                for(auto it = loop.begin()+2; it != loop.end(); ++it) {
                    S.col(1) = S.col(2);
                    S.col(2) = V.col(*it);
                    vol += mtao::geometry::volume_signed(S);
                }
                return vol;
                //return mtao::geometry::volume_signed(S);
            };
            for(auto&& loop: indices) {
                vol += loop_volume(loop);
            }
            return vol;
        }
}
