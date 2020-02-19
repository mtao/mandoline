#include "mandoline/mesh2.hpp"
#include "mandoline/construction/generator.hpp"

namespace mandoline::construction {
    template <>
        class CutCellGenerator<2>: public CutCellEdgeGenerator<2> {
            public:
                using CCEG = CutCellEdgeGenerator<2>;
                using BoundaryElements = typename CCEG::BoundaryElements;
                using CCEG::CCEG;
                using CoordType = typename CCEG::CoordType;
                CutCellGenerator() = default;
                CutCellGenerator(CutCellGenerator&&) = default;
                CutCellGenerator& operator=(CutCellGenerator&&) = default;
                CutCellMesh<2> generate() const ;
                void add_boundary_elements(const BoundaryElements& E);
                void bake_faces() override;
                void extra_metadata(CutCellMesh<2>& mesh) const;
                mtao::geometry::mesh::HalfEdgeMesh hem;

                mtao::map<int,CutFace<D>> m_faces;
                std::set<int> mesh_face_indices;
                std::array<std::set<Edge>,3> axial_primal_faces;
                mtao::map<int,int> cut_cell_to_primal_map; // store the cut-face -> input face map
        };


    template <>
        CutCellMesh<2> CutCellEdgeGenerator<2>::generate_faces() const;
    template <>
        CutCellMesh<2> CutCellEdgeGenerator<2>::generate() const;
    template <>
        auto CutCellEdgeGenerator<2>::compute_planar_hem(const ColVecs& V, const Edges& E, const GridDatab& interior_cell_mask, int cell_size) const-> std::tuple<mtao::geometry::mesh::HalfEdgeMesh,std::set<Edge>>;
    template <>
        auto CutCellEdgeGenerator<2>::compute_planar_hem(const std::vector<VType>& V, const Edges& E, const GridDatab& interior_cell_mask, int cell_size) const-> std::tuple<mtao::geometry::mesh::HalfEdgeMesh,std::set<Edge>>;
    template <>
        auto CutCellEdgeGenerator<2>::compute_planar_hem(const std::vector<VType>& GV, const ColVecs& V, const Edges& E, const GridDatab& interior_cell_mask, int cell_size) const-> std::tuple<mtao::geometry::mesh::HalfEdgeMesh,std::set<Edge>>;
}
