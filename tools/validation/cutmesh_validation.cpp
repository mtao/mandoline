#include <iostream>
#include "mandoline/mesh3.hpp"
#include <mtao/cmdline_parser.hpp>
#include <mtao/logging/logger.hpp>
#include "cutmesh_validation.hpp"
using namespace mtao::logging;

using namespace mandoline;


int main(int argc, char * argv[]) {

    auto&& log = make_logger("profiler",mtao::logging::Level::All);
    mtao::CommandLineParser clp;
    clp.parse(argc, argv);

    if(clp.args().size() < 1) {
        fatal() << "No input mesh filename!";
        return {};
    }

    std::string input_cutmesh = clp.arg(0);

    CutCellMesh<3> ccm = CutCellMesh<3>::from_proto(input_cutmesh);

    ccm.triangulate_faces();
    bool pcwn = pcwn_check(ccm);
    bool is_paired = paired_boundary(ccm);
    bool ffu = faces_fully_utilized(ccm);
    bool grid_cells = grid_cells_fully_utilized(ccm);
    bool region_counts = equal_region_check(ccm);
    bool volume_valid  = volume_check(ccm);

    /*
    "name","pcwn","is_paired","ffu","grid_cell","region_counts","volume_valid" 
    */
    mtao::logging::fatal() << clp.arg(0) << ","  <<
    pcwn << "," << 
    is_paired << "," << 
    ffu << "," << 
    grid_cells << "," << 
    region_counts << "," << 
    volume_valid ;

    //mtao::logging::fatal() << clp.arg(0) << "," << (ccm.folded_faces.size() == 0) << "," << (is_not_pcwn == 0) << "," << is_paired << "," << ffu;
    //mtao::logging::fatal() << clp.arg(0) << "," << (ccm.folded_faces.size() == 0) << "," << is_pcwn << "," << is_not_pcwn << "," << is_paired << "," << ffu;




    return 0;
}


