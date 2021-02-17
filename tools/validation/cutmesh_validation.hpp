#include "mandoline/mesh3.hpp"

bool run_all(const mandoline::CutCellMesh<3>& ccm);
// grid cell isolation
// pcwn
// number of cells taht are pcwn {{is_pcwn,is_not_pcwn}}
std::array<int, 2> pcwn_count(const mandoline::CutCellMesh<3>&);
bool pcwn_check(const mandoline::CutCellMesh<3>&);
// surface area
bool faces_fully_utilized(const mandoline::CutCellMesh<3>&);
// grid volume
// number of grid cells underutilized and number of grid cells completely
// mismatched (either empty or used despite not being masked)
std::array<int, 2> grid_cells_fully_utilized_count(
    const mandoline::CutCellMesh<3>& ccm);
bool grid_cells_fully_utilized(const mandoline::CutCellMesh<3>&);
// region counts {{mandoline regions, igl regions}}
std::array<int, 2> region_counts(const mandoline::CutCellMesh<3>&);
bool equal_region_check(const mandoline::CutCellMesh<3>&);
// input volume
bool volume_check(const mandoline::CutCellMesh<3>&);

// face topological utilization
bool paired_boundary(const mandoline::CutCellMesh<3>&);

// exterior cells should have valence 6
bool exterior_cell_valence_counts(const mandoline::CutCellMesh<3>&);

// UTILITY
//
mtao::ColVecs2i regions(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F);
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> orig_with_bbox(
    const mandoline::CutCellMesh<3>&);
std::array<int, 2> region_counts(const mandoline::CutCellMesh<3>&,
                                 const mtao::ColVecs2i&);
mtao::ColVecs2i input_mesh_regions(const mandoline::CutCellMesh<3>&);
std::map<int, double> region_volumes(const mandoline::CutCellMesh<3>&);
mtao::VecXd brep_region_volumes(const mandoline::CutCellMesh<3>&);
Eigen::VectorXd brep_region_volumes(const mtao::ColVecs3d& V,
                                    const mtao::ColVecs3i& F);
Eigen::VectorXd brep_region_volumes(const mtao::ColVecs3d& V,
                                    const mtao::ColVecs3i& F,
                                    const mtao::ColVecs2i& C);
// for a set of vertices that comprise a face, list the set of potential cells
// that it resides within
std::set<std::array<int, 3>> possible_cells(
    const mandoline::CutCellMesh<3>& ccm, const std::vector<int>& face);

// each exterior grid entry should have valence 6
std::map<int, int> exterior_grid_valences(const mandoline::CutCellMesh<3>&);

bool cut_normals_match_parents(const mandoline::CutCellMesh<3>& ccm);
// bool cut_tangents_match_parents(const mandoline::CutCellMesh<2>& ccm);
bool cut_tangents_match_parents(const mandoline::CutCellMesh<3>& ccm);
