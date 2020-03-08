#include "mandoline/construction/generator3.hpp"
#include <cxxopts.hpp>


mandoline::construction::CutCellGenerator<3> make_generator(const cxxopts::ParseResult& clp);
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> read_mesh_input(const cxxopts::ParseResult& clp);
mandoline::construction::CutCellGenerator<3> make_generator(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const cxxopts::ParseResult& clp);
mandoline::CutCellMesh<3> make_cutmesh(const mandoline::construction::CutCellGenerator<3>& ccg, const cxxopts::ParseResult& clp);
