#include "mandoline/construction/generator.hpp"
#include <mtao/cmdline_parser.hpp>


mandoline::construction::CutCellGenerator<3> make_generator(const mtao::CommandLineParser& clp);
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> read_mesh_input(const mtao::CommandLineParser& clp);
mandoline::construction::CutCellGenerator<3> make_generator(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const mtao::CommandLineParser& clp);
mandoline::CutCellMesh<3> make_cutmesh(const mandoline::construction::CutCellGenerator<3>& ccg, const mtao::CommandLineParser& clp);
