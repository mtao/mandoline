#include "mandoline/mesh3.hpp"
#include <mtao/cmdline_parser.hpp>
#include <cxxopts.hpp>


cxxopts::Options make_cutmesh_clparser();
mandoline::CutCellMesh<3> make_cutmesh(const cxxopts::ParseResult& clp);
mandoline::CutCellMesh<3> make_cutmesh(int argc, char * argv[]);
