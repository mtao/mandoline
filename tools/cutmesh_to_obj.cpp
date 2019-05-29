#include "mandoline/mesh3.hpp"
#include <mtao/cmdline_parser.hpp>
#include <mtao/logging/logger.hpp>
using namespace mtao::logging;
using namespace mandoline;



int main(int argc, char * argv[]) {

    auto&& log = make_logger("profiler",mtao::logging::Level::All);
    mtao::CommandLineParser clp;
    clp.add_option("open-regions",false);
    clp.add_option("write-regions",true);
    clp.add_option("write-monolithic",true);
    clp.add_option("write-separate",true);
    clp.add_option("normalize",false);
    clp.add_option("normalize-unit",false);
    clp.add_option("cell-grid-ownership",false);
    clp.parse(argc, argv);

    if(clp.args().size() < 1) {
        fatal() << "No input mesh filename!";
        return {};
    }
    if(clp.args().size() < 2) {
        fatal() << "No output mesh filename!";
    }

    std::string input_cutmesh = clp.arg(0);
    std::string output_prefix = clp.arg(1);
    bool write_regions = clp.optT<bool>("write-regions");
    bool normalize = clp.optT<bool>("normalize");
    bool normalize_gc = clp.optT<bool>("normalize-unit");
    bool write_monolithic = clp.optT<bool>("write-monolithic");
    bool write_separate = clp.optT<bool>("write-separate");
    bool open_regions = clp.optT<bool>("open-regions");
    bool cell_grid_ownership = clp.optT<bool>("cell-grid-ownership");

    CutCellMesh<3> ccm = CutCellMesh<3>::from_proto(input_cutmesh);




    ccm.triangulate_faces();


    if(write_regions   ) {  
        ccm.write_obj_regions(output_prefix);
    }
    if(write_monolithic) {  
        ccm.write_obj(output_prefix);
    }
    if(write_separate  ) { 
        ccm.write_obj_separate(output_prefix,false);
    }
    std::cout << "woop" << std::endl;
    if(cell_grid_ownership) {
        std::cout << "Woop" << std::endl;

        std::ofstream ofs(output_prefix + "_cells.txt");
        auto C = ccm.cell_centroids();
        for(auto&& [i,cell]: mtao::iterator::enumerate(ccm.cells())) {
            ofs << i << " ";
            auto cent = C.col(i);
            auto [c,q] = ccm.coord(cent);
            std::copy(c.begin(),c.end(),std::ostream_iterator<int>(ofs," "));
            ofs << std::endl;
        }
        for(auto&& [i,c]: ccm.adaptive_grid().cells) {
            auto center  = c.corner();
            int w = c.width();
            for(auto&& c: center) {
                c += w/2;
            }
            ofs << i << " ";
            std::copy(center.begin(),center.end(),std::ostream_iterator<int>(ofs," "));
            ofs << std::endl;
        }
    }
    return 0;
}


