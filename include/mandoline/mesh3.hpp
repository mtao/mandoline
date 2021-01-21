#pragma once
#include "mandoline/barycentric_triangle_face.hpp"
#include "mandoline/cutcell.hpp"
#include "mandoline/cutface.hpp"
#include "mesh.hpp"
#if defined(MANDOLINE_USE_ADAPTIVE_GRID)
#include "mandoline/adaptive_grid.hpp"
#else
#include "mandoline/exterior_grid.hpp"
#endif
namespace mandoline {
namespace construction {
template <int D>
class CutCellGenerator;
template <>
class CutCellGenerator<3>;
}  // namespace construction
template <>
struct CutCellMesh<3> : public CutCellMeshBase<3, CutCellMesh<3>> {
    // NOTE: Grid index of -1 == inside stencil, -2 == boundary
   public:
    friend class construction::CutCellGenerator<3>;
    using Base = CutCellMeshBase<3, CutCellMesh<3>>;
    using coord_type = typename Base::coord_type;
    using ColVecs = typename Base::ColVecs;
    using VecX = typename Base::VecX;
    using Faces = mtao::ColVectors<int, 3>;
    using Face = mtao::Vector<int, 3>;
    using Vertex = typename Base::Vertex;
#if defined(MANDOLINE_USE_ADAPTIVE_GRID)
    using ExteriorGridType = AdaptiveGrid;
#else
    using ExteriorGridType = ExteriorGrid<3>;
#endif

    CutCellMesh() = default;
    CutCellMesh(const CutCellMesh &) = default;
    CutCellMesh(CutCellMesh &&) = default;
    CutCellMesh &operator=(const CutCellMesh &) = default;
    CutCellMesh &operator=(CutCellMesh &&) = default;
    CutCellMesh(const StaggeredGrid &g, const std::vector<Vertex> &v = {})
        : Base(g, v), m_exterior_grid(g) {}

    static CutCellMesh<3> from_file(const std::string &filename);

    // Inherited edges() function doesn't take advantage of the staggered grid
    Edges edges() const;
    const std::vector<CutFace<3>> &faces() const { return m_faces; }
    const std::vector<CutFace<3>> &cut_faces() const { return m_faces; }
    const CutFace<3> &cut_face(size_t index) const { return m_faces.at(index); }
    const auto &cells() const { return m_cells; }
    const ExteriorGridType &exterior_grid() const { return m_exterior_grid; }
    const mtao::ColVecs3i &origF() const { return m_origF; }
    const mtao::map<int, BarycentricTriangleFace> &mesh_cut_faces() const {
        return m_mesh_cut_faces;
    }

    // info on cells
    size_t cell_size() const;
    size_t num_cells() const;
    size_t cut_cell_size() const;
    size_t num_cut_cells() const;
    bool is_cut_cell(int index) const;
    bool is_exterior_cell(int index) const;
    std::vector<int> regions(bool boundary_sign_regions = false) const;
    std::vector<std::array<std::set<int>, 2>> face_regions() const;
    std::vector<std::array<std::set<int>, 2>> orig_face_regions() const;
    mtao::ColVecs3d region_centroids() const;
    std::map<coord_type, std::set<int>> cells_by_grid_cell() const;
    std::set<int> cells_in_grid_cell(const coord_type &c) const;
    int get_cell_index(const VecCRef &p) const;

    bool is_in_cell(const VecCRef &p, size_t index) const;

    // info on faces
    size_t face_size() const;
    size_t cut_face_size() const;
    size_t num_faces() const;
    size_t num_cut_faces() const;
    bool is_cut_face(size_t index) const;
    bool is_mesh_face(int idx) const;
    std::vector<bool> boundary_faces() const;

    // grid cell mask / gneeration info
    const GridDatab &active_cell_mask() const;
    std::set<coord_type> active_cells() const;
    size_t active_cell_count() const;
    // grid cell index from a local coordinate (grid is on integer coordinates)
    int local_grid_cell_index(const VecCRef &) const;
    // grid cell index from a world space coordinate
    int world_grid_cell_index(const VecCRef &) const;

    // Differential geometry info
    // cell -> face boundary operator
    Eigen::SparseMatrix<double> boundary(
        bool include_domain_boundary_faces = false) const;
    mtao::ColVecs3d face_centroids() const;
    mtao::ColVecs3d cell_centroids() const;
    ColVecs dual_vertices() const;  // alias for centroids
    VecX cell_volumes() const;
    VecX face_volumes(bool from_triangulation = false) const;

    // impl in operators/volume.hpp
    mtao::VecXd dual_edge_lengths() const;
    mtao::VecXd dual_hodge2() const;
    mtao::VecXd primal_hodge2() const;
    mtao::VecXd dual_hodge3() const;
    mtao::VecXd primal_hodge3() const;
    // impl in operators/masks.hpp
    mtao::VecXd mesh_face_mask() const;      // for removing mesh faces
    mtao::VecXd grid_boundary_mask() const;  // for removing mesh faces

    std::set<int> grid_boundary_faces() const;  // for removing mesh faces

    // Eigen::SparseMatrix<double> trilinear_matrix() const;

    // impl in operators/boundary.hpp
    // mesh vertex -> cut vertex
    Eigen::SparseMatrix<double> barycentric_matrix() const;
    // mesh face -> cut face
    Eigen::SparseMatrix<double> face_barycentric_volume_matrix() const;
    // grid vertex -> cut vertex
    Eigen::SparseMatrix<double> trilinear_matrix() const;
    // grid face -> cut face
    Eigen::SparseMatrix<double> face_grid_volume_matrix() const;

    // an unsigned boundary map when its not necessary
    std::vector<std::set<int>> cell_faces() const;
    std::set<int> cell_faces(int idx) const;

    // serialization
    void write(const std::string &filename) const;
    void serialize(protobuf::CutMeshProto &) const;
    static CutCellMesh<3> from_proto(const protobuf::CutMeshProto &);
    static CutCellMesh<3> from_proto(const std::string &filename);

    // Caches triangulations for each CutFace, important for triangulating
    // things like cells
    void triangulate_faces(bool add_verts = true);

    // If the input ColVecs3d has nonzero size then the mesh is with reference
    // to those vertices Triangulation of different mesh elements
    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> triangulate_face(
        int face_index) const;
    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> triangulated_cell(
        int cell_index, bool use_base = true, bool use_flap = true) const;

    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> compact_triangulated_cell(
        int cell_index) const;
    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> compact_triangulated_face(
        int face_index, bool invert = false) const;

    // Vertices as projected into different 2d planes
    std::array<mtao::ColVecs2d, 3> compute_subVs() const;

#if defined(MANDOLINE_USE_ADAPTIVE_GRID)
    const std::map<int, int> &adaptive_grid_regions() const {
        return m_adaptive_grid_regions;
    }
#endif
    const std::set<int> &folded_faces() const { return m_folded_faces; }
    bool is_folded_face(int idx) const {
        return m_folded_faces.find(idx) != m_folded_faces.end();
    }

   private:
    // Primary geometry data
    std::vector<CutFace<3>> m_faces;
    std::vector<CutCell> m_cells;

    ExteriorGridType m_exterior_grid;
#if defined(MANDOLINE_USE_ADAPTIVE_GRID)
    // Cell annotations
    std::map<int, int> m_adaptive_grid_regions;
#endif

    // a map from cut-faces to their intrinsic representation on a mesh face
    mtao::map<int, BarycentricTriangleFace> m_mesh_cut_faces;

    // Original mesh
    mtao::ColVecs2i m_origE;
    mtao::ColVecs3i m_origF;

    // Face annotations
    std::array<std::set<int>, 3> m_axial_faces;
    std::set<int> m_folded_faces;

   private:
    void recompute_active_cell_mask();
    void write_obj(const std::string &prefix, const std::set<int> &indices,
                   const std::optional<int> &region = {},
                   bool show_indices = false, bool show_base = true,
                   bool show_flaps = false, bool mesh_face = false) const;
    // mtao::ColVecs3i origF;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace mandoline
