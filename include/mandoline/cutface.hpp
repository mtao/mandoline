#pragma once
#include "mandoline/coord_masked_geometry.hpp"
#include "mandoline/vertex.hpp"
#include <mtao/geometry/grid/grid.h>
#include <igl/solid_angle.h>

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
            CutFace(const coord_mask<D>& m, const std::set<std::vector<int>>& indices, const std::array<int,2>& id): Base(m,indices), id(id),N(mtao::Vec3d::Unit(m.bound_axis())) {}

            CutFace(const CutMeshFace<D>& F, const mtao::Vec3d& N): CutFace(F.mask(),F.indices,F.parent_fid,N) {}
            CutFace() = default;
            CutFace(const CutFace&) = default;
            CutFace(CutFace&&) = default;
            CutFace& operator=(const CutFace&) = default;
            CutFace& operator=(CutFace&&) = default;

            //a tuple where the first entry is the cell-grid index of hte external cell, and the bool is whether it is below it or not 
            //  i.e given two cells cm cp bounded by f, 
            //  [cm] f [cp], if [cm] is not part of the staggered grid then we store
            //  cm,1
            //  otherwise we store 
            //  cp,0
            // this sign is with respect to "the boundary operator of hte cell is of positive sign", which may be backwards
            std::optional<std::tuple<int,bool>> external_boundary = {};
            std::optional<mtao::ColVecs3d> triangulated_vertices;
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

        template <typename Derived, typename VecType>
            double solid_angle(const Eigen::MatrixBase<Derived>& V, const Eigen::MatrixBase<VecType>& v) const;

            void  serialize(protobuf::CutFace& face) const ;
            static CutFace<D> from_proto(const protobuf::CutFace& face) ;
            template <typename Derived>
            void update_mask(const std::vector<Vertex<D>>& V, const mtao::geometry::grid::indexing::IndexerBase<3,Derived>& indexer);

            mtao::ColVecs3i triangulate_fan() const;
            mtao::ColVecs3i triangulate_earclipping(const mtao::ColVecs2d& V) const;
            std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> triangulate_triangle(const mtao::ColVecs2d& V, bool add_vertices = false) const;
            mtao::ColVecs3i triangulate(const std::array<mtao::ColVecs2d,D>& V) const;
            std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> triangulate(const std::array<mtao::ColVecs2d,D>& V, bool add_vertices) const;
            std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> triangulate(const mtao::ColVecs2d& V, bool add_vertices) const;
            void cache_triangulation(const std::array<mtao::ColVecs2d,D>& V, bool add_verts = true) ;
            void cache_triangulation(const mtao::ColVecs3i& F) ;
            void cache_triangulation(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F) ;







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
        std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> CutFace<3>::triangulate_triangle(const mtao::ColVecs2d& V, bool add_vertices) const;
    template <>
        mtao::ColVecs3i CutFace<3>::triangulate(const std::array<mtao::ColVecs2d,3>& V) const;
    template <>
        std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> CutFace<3>::triangulate(const std::array<mtao::ColVecs2d,3>& V, bool) const;
    template <>
        std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> CutFace<3>::triangulate(const mtao::ColVecs2d& V, bool) const;
    template <>
        void CutFace<3>::cache_triangulation(const std::array<mtao::ColVecs2d,3>& V, bool add_verts) ;
    template <>
        void CutFace<3>::cache_triangulation(const mtao::ColVecs3i& F) ;
    template <>
        void CutFace<3>::cache_triangulation(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F) ;

    template <>
        mtao::ColVecs3i CutMeshFace<3>::triangulate() const;

    template <>
        template <typename Derived, typename VecType>
            double CutFace<3>::solid_angle(const Eigen::MatrixBase<Derived>& V, const Eigen::MatrixBase<VecType>& v) const {

                double sa = 0;
                if(triangulation) {
                    auto& F = *triangulation;
                    if(triangulated_vertices) {
                        auto& V = *triangulated_vertices;
                        for(int i = 0; i < F.cols(); ++i) {
                            auto f = F.col(i);

                            sa += igl::solid_angle(
                                    V.col(f(0)),V.col(f(1)),V.col(f(2))
                                    ,v);

                        }
                    } else {
                        for(int i = 0; i < F.cols(); ++i) {
                            auto f = F.col(i);

                            sa += igl::solid_angle(
                                    V.col(f(0)),V.col(f(1)),V.col(f(2))
                                    ,v);

                        }
                    }
                } else {
                    for(auto&& f: indices) {
                        //just use triangle fan!
                        for(int i = 1; i < f.size()-1; ++i)
                        {
                            sa += igl::solid_angle(V.col(f[0]),V.col(f[i]),V.col(f[i+1]),v);
                        }
                    }
                }
                return sa;
            }

}


#include "mandoline/cutface_impl.hpp"
