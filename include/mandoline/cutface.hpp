#pragma once
#include "mandoline/coord_masked_geometry.hpp"

namespace mandoline {

    template <int D>
        struct CutMeshFace: public CoordMaskedSimplePolygon<D> {

            using Base = CoordMaskedSimplePolygon<D>;
            using Base::mask;
            using Base::indices;
            using Base::Base;


            CutMeshFace(const coord_mask<D>& m, const std::vector<int>& indices, int pid): Base(m,indices), parent_fid(pid) {}
            CutMeshFace() = default;
            CutMeshFace(const CutMeshFace&) = default;
            CutMeshFace(CutMeshFace&&) = default;
            CutMeshFace& operator=(const CutMeshFace&) = default;
            CutMeshFace& operator=(CutMeshFace&&) = default;
            size_t size() const { return indices.size(); }
            mtao::ColVecs3i triangulate() const;

            int parent_fid = -1;

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        };

    template <int D>
        struct CutFace: public CoordMaskedPolygon<D> {

            using Base = CoordMaskedPolygon<D>;
            using Base::mask;
            using Base::indices;
            using Base::Base;

            using IDType = std::variant<int,std::array<int,2>>;


            CutFace(const coord_mask<D>& m, const std::vector<int>& indices, const IDType& id, const mtao::Vec3d& N): Base(m,{indices}), id(id), N(N) {}
            CutFace(const coord_mask<D>& m, const std::set<std::vector<int>>& indices, const std::array<int,2>& id): Base(m,indices), id(id),N(mtao::Vec3d::Unit(m.unbound_axis())) {}

            CutFace(const CutMeshFace<D>& F, const mtao::Vec3d& N): CutFace(F.mask(),F.indices,F.parent_fid,N) {}
            CutFace() = default;
            CutFace(const CutFace&) = default;
            CutFace(CutFace&&) = default;
            CutFace& operator=(const CutFace&) = default;
            CutFace& operator=(CutFace&&) = default;

            std::optional<std::tuple<int,bool>> external_boundary = {};
            std::optional<mtao::ColVecs3i> triangulation;
            bool is_mesh_face() const ;
            bool is_axial_face() const ;
            int as_face_id() const ;
            const std::array<int,2>& as_axial_id() const ;
            int as_axial_axis() const ;
            int as_axial_coord() const ;


            operator std::string() const ;

            template <typename Derived>
                mtao::Vec3d brep_centroid(const Eigen::MatrixBase<Derived>& V, bool use_triangulation = false) const ;

            template <typename Derived>
                double brep_volume(const Eigen::MatrixBase<Derived>& V, bool use_triangulation = false) const ;


            void  serialize(CutMeshProto::CutFace& face) const ;
            static CutFace<D> from_proto(const CutMeshProto::CutFace& face) ;

            mtao::ColVecs3i triangulate_fan() const;
            mtao::ColVecs3i triangulate_earclipping(const mtao::ColVecs2d& V) const;
            mtao::ColVecs3i triangulate_triangle(const mtao::ColVecs2d& V) const;
            mtao::ColVecs3i triangulate(const std::array<mtao::ColVecs2d,D>& V) const;
            void cache_triangulation(const std::array<mtao::ColVecs2d,D>& V) ;
            void cache_triangulation(const mtao::ColVecs3i& F) ;







            IDType id;
            mtao::Vec3d N;
            public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        };


    template <>
        mtao::ColVecs3i CutFace<3>::triangulate_fan() const;
    template <>
        mtao::ColVecs3i CutFace<3>::triangulate_earclipping(const mtao::ColVecs2d& V) const;
    template <>
        mtao::ColVecs3i CutFace<3>::triangulate_triangle(const mtao::ColVecs2d& V) const;
    template <>
        mtao::ColVecs3i CutFace<3>::triangulate(const std::array<mtao::ColVecs2d,3>& V) const;
    template <>
        void CutFace<3>::cache_triangulation(const std::array<mtao::ColVecs2d,3>& V) ;
    template <>
        void CutFace<3>::cache_triangulation(const mtao::ColVecs3i& F) ;

    template <>
        mtao::ColVecs3i CutMeshFace<3>::triangulate() const;

}


#include "mandoline/cutface_impl.hpp"
