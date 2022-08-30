#include "mandoline/mesh2.hpp"

namespace mandoline::tools {

// Input:
// Vertices
// vertex_idx, xcoord, ycoord
// vertex_idx, xcoord, ycoord
//
// Edges
// edge_idx, vertex_idx1, vertex_idx2, left_element_idx, right_element_idx
// edge_idx, vertex_idx1, vertex_idx2, left_element_idx, right_element_idx
std::tuple<balsa::eigen::ColVecs2d, balsa::eigen::ColVecs4i> read_plcurve(const std::string &filename);


void write_plcurve(const balsa::eigen::ColVecs2d &V, const balsa::eigen::ColVecs2i &E, const std::string &filename);
void write_plcurve(const balsa::eigen::ColVecs2d &V, const balsa::eigen::ColVecs4i &E, const std::string &filename);


void write_cutmesh2_plcurve(const mandoline::CutCellMesh<2> &ccm, const std::string &filename, const balsa::eigen::ColVecs4i &inputE = {});
}// namespace mandoline::tools
