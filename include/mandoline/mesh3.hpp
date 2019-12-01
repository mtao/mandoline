#pragma once
#include "mesh.hpp"
namespace mandoline {
    namespace construction {
        template <int D>
        class CutCellGenerator;
        template <>
        class CutCellGenerator<3>;
    }
    template <>
        struct CutCellMesh<3>: public CutCellMeshBase<3,CutCellMesh<3>> {
            public:
                friend class construction::CutCellGenerator<3>;
                using Base = CutCellMeshBase<3,CutCellMesh<3>>;
                using coord_type = typename Base::coord_type;
                using ColVecs = typename Base::ColVecs;
                using VecX = typename Base::VecX;
                using Faces = mtao::ColVectors<int,3>;
                using Face = mtao::Vector<int,3>;
                using Vertex = typename Base::Vertex;

                CutCellMesh() = default;
                CutCellMesh(const CutCellMesh&) = default;
                CutCellMesh(CutCellMesh&&) = default;
                CutCellMesh& operator=(const CutCellMesh&) = default;
                CutCellMesh& operator=(CutCellMesh&&) = default;
                CutCellMesh(const StaggeredGrid& g, const std::vector<Vertex>& v = {}): Base(g,v), m_adaptive_grid(g) {}

                static CutCellMesh<3> from_file(const std::string& filename);

                const auto& faces() const { return m_faces; }
                const auto& cells() const { return m_cells; }
                const AdaptiveGrid& adaptive_grid() const { return m_adaptive_grid; }
                const mtao::ColVecs3d& origV() const { return m_origV; }
                const mtao::ColVecs3i& origF() const { return m_origF; }


                //info on cells
                size_t cell_size() const;
                size_t cut_cell_size() const;
                bool is_cut_cell(int index) const;
                bool is_adaptive_cell(int index) const ;
                std::vector<int> regions(bool boundary_sign_regions=false) const;
                std::vector<std::array<std::set<int>,2>> face_regions() const;
                std::vector<std::array<std::set<int>,2>> orig_face_regions() const;
                mtao::ColVecs3d region_centroids() const;
                std::map<coord_type,std::set<int>> cells_by_grid_cell() const;
                std::set<int> cells_in_grid_cell(const coord_type& c) const;
                int get_cell_index(const VecCRef& p) const;

                //info on faces
                size_t face_size() const;
                size_t cut_face_size() const;
                bool is_mesh_face(int idx) const;
                std::vector<bool> boundary_faces() const;

                //grid cell mask / gneeration info
                GridDatab active_cell_mask() const;
                std::set<coord_type> active_cells() const;
                size_t active_cell_count() const;
                //grid cell index from a local coordinate (grid is on integer coordinates)
                int local_grid_cell_index(const VecCRef&) const;
                //grid cell index from a world space coordinate
                int world_grid_cell_index(const VecCRef&) const;

                //Differential geometry info
                Eigen::SparseMatrix<double> boundary() const;
                VecX cell_volumes() const;
                mtao::ColVecs3d face_centroids() const;
                mtao::ColVecs3d cell_centroids() const;
                ColVecs dual_vertices() const;//alias for centroids
                VecX face_volumes(bool from_triangulation = false)const ;
                mtao::VecXd dual_edge_lengths() const;
                mtao::VecXd dual_hodge2() const;
                mtao::VecXd primal_hodge2() const;
                mtao::VecXd dual_hodge3() const;
                mtao::VecXd primal_hodge3() const;
                mtao::VecXd mesh_face_mask() const;//for removing mesh faces

                //Eigen::SparseMatrix<double> trilinear_matrix() const;

                //mesh vertex -> cut vertex
                Eigen::SparseMatrix<double> barycentric_matrix() const;
                //mesh face -> cut face
                Eigen::SparseMatrix<double> face_barycentric_volume_matrix() const;
                //grid vertex -> cut vertex
                Eigen::SparseMatrix<double> trilinear_matrix() const;
                //grid face -> cut face
                Eigen::SparseMatrix<double> face_grid_volume_matrix() const;

                //an unsigned boundary map when its not necessary
                std::vector<std::set<int>> cell_faces() const;
                std::set<int> cell_faces(int idx) const;

                //serialization
                void write(const std::string& filename) const;
                void  serialize(CutMeshProto&) const;
                static CutCellMesh<3> from_proto(const CutMeshProto&);
                static CutCellMesh<3> from_proto(const std::string& filename);


                //Caches triangulations for each CutFace, important for triangulating things like cells
                void triangulate_faces();

                //If the input ColVecs3d has nonzero size then the mesh is with reference to those vertices
                //Triangulation of different mesh elements
                std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> triangulate_face(int face_index) const;
                std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> triangulated_cell(int cell_index, bool use_base = true, bool use_flap = true) const;

                std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> compact_triangulated_cell(int cell_index) const;
                std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> compact_triangulated_face(int face_index) const;



                //Triangle mesh writers
                void write_obj(const std::string& filename) const;
                void write_obj_flaps(const std::string& filename) const;
                void write_obj_regions(const std::string& filename) const;
                void write_obj_separate(const std::string& prefix, bool separate_flaps = false) const;
                void write_mesh_obj_separate(const std::string& prefix) const;
                void write_mesh_surface_obj(const std::string& prefix) const;







                //Vertices as projected into different 2d planes
                std::array<mtao::ColVecs2d,3> compute_subVs() const;



                const std::map<int,int>& adaptive_grid_regions() const { return m_adaptive_grid_regions; }
                const std::set<int>& folded_faces() const { return m_folded_faces;}
                bool is_folded_face(int idx) const { return m_folded_faces.find(idx) != m_folded_faces.end();}



            private:
                //Primary geometry data
                std::vector<CutFace<3>> m_faces;
                std::vector<CutCell> m_cells;
                AdaptiveGrid m_adaptive_grid;

                //support for mapping
                mtao::map<int,BarycentricTriangleFace> m_mesh_faces;

                //Original mesh
                mtao::ColVecs3d m_origV;
                mtao::ColVecs3i m_origF;

                //Face annotations
                std::array<std::set<int>,3> m_axial_faces;
                std::set<int> m_folded_faces;

                //Cell annotations
                std::map<int,int> m_adaptive_grid_regions;


            private:
                void write_obj(const std::string& prefix, const std::set<int>& indices, const std::optional<int>& region = {}, bool show_indices = false, bool show_base=true, bool show_flaps = false, bool mesh_face = false ) const;
                //mtao::ColVecs3i origF;
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        };


}
