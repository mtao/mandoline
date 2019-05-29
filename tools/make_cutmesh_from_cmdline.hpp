#include "mandoline/mesh3.hpp"
#include <mtao/cmdline_parser.hpp>


mandoline::CutCellMesh<3> make_cutmesh(const mtao::CommandLineParser& clp);
mtao::CommandLineParser make_cutmesh_clparser();
mandoline::CutCellMesh<3> make_cutmesh(int argc, char * argv[]);
