#include "mandoline/construction/generator.hpp"
#include "mandoline/mesh3.hpp"
#include "mandoline/adaptive_grid.hpp"

namespace mandoline::construction {
template<>
class CutCellGenerator<3> : public CutCellEdgeGenerator<3> {
  public:
    using CCEG = CutCellEdgeGenerator<3>;
    using BoundaryElements = typename CCEG::BoundaryElements;
    using coord_type = typename CCEG::coord_type;
    using CCEG::add_boundary_elements;
    using CCEG::CCEG;
    CutCellGenerator() = default;
    CutCellGenerator(CutCellGenerator &&) = default;
    CutCellGenerator &operator=(CutCellGenerator &&) = default;
    using CCEG::vertex_shape;

    size_t Vsize() const {
        return CCEG::Vsize();
    }
    ~CutCellGenerator();

    const mtao::map<int, CutFace<D>> &faces() const { return m_faces; }

    void compute_faces();
    void compute_faces_vertex();
    void compute_faces_axis(int idx);
    mtao::map<int, CutFace<D>> compute_faces_axis(int idx, int cidx) const;

    void bake_faces() override;
    void bake_cells() override;
    void update_active_grid_cell_mask();
#if defined(MANDOLINE_USE_ADAPTIVE_GRID)
    std::optional<int> adaptive_level = 0;
    std::optional<AdaptiveGrid> adaptive_grid;
    std::optional<std::map<int, int>> adaptive_grid_regions;
#else
    std::optional<ExteriorGrid<3>> exterior_grid;;
#endif
    void set_region_map(const std::map<int, std::set<int>> &region_vertices);
    void bake() override;
    void clear() override;


    static mtao::ColVecs3i faceMap_to_faces(mtao::map<int, mtao::ColVecs3i> &fm);

    void update_vertices_from_intersections();
    CutCellMesh<3> generate() const;
    //void make_faces(const CutCellMesh<3>& ccm ) const;

    void extra_metadata(CutCellMesh<3> &mesh) const;
    /*
                   const auto& newV() const { return m_newV; }
                   const Vec& newV(int i) const { return m_newV[i]; }
                   auto&& stl_face_cutE() const { return m_face_cutedges; }
                   const BoundaryElements& origF() const { return m_origF;}
                   auto origF(int i) const { return m_origF.col(i);}
                   const BoundaryElements& origFE() const { return m_cell_edge_map;}
                   auto origFE(int i) const { return m_cell_edge_map.col(i);}
                   */

    //baked by bake_faces
    struct AxisHEMData {
        using GridDatab = mtao::geometry::grid::GridDataD<bool, 2>;
        GridDatab active_grid_cell_mask;
        std::set<Edge> edges;
        std::set<Edge> boundary_edges;
        mtao::geometry::mesh::HalfEdgeMesh hem;
        bool is_boundary_cell(const std::vector<int> &verts) const;
    };
    bool check_cell_containment() const;
    bool check_face_utilization() const;

    mtao::Vec3d area_normal(const std::vector<int> &F) const;
    mtao::Vec3d area_normal(const std::set<std::vector<int>> &F) const;
    std::set<Edge> edge_slice(int dim, int slice) const;
    mtao::map<int, int> cut_cell_to_primal_map;// store the cut-face -> input face map
    mtao::ColVecs3d origN;// the normals from the input mesh
    std::array<mtao::map<int, AxisHEMData>, 3> axis_hem_data;// per-cut-plane information
    std::array<std::set<Edge>, 3> axial_primal_faces;// a hash for the faces of the input mesh that lie on axial planes
    BoundaryElements m_newF;// ??? what is this

    mtao::map<int, CutFace<D>> m_faces;
    std::set<int> mesh_face_indices;// as m_faces loses track of teh cutmesh faces, this keeps track
    std::array<std::set<int>, 3> axis_face_indices;// store the faces that come from each axis
    std::set<int> folded_faces;// elements that are on the boundary of hte input mesh


    std::vector<CutCell> cell_boundaries;
    std::set<int> boundary_vertices;
    std::set<int> boundary_faces;
};

template<>
CutCellMesh<3> CutCellEdgeGenerator<3>::generate() const;
}// namespace mandoline::construction
