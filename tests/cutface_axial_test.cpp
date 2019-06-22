#include <mandoline/mesh3.hpp>
#include <mtao/cmdline_parser.hpp>
#include <mtao/logging/timer.hpp>
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


    for(auto&& face: ccm.faces()) {
        if(face.is_axial_face()) {
            if(face.count() != 1) {
                std::cout << "Bad mask on face!" << ": " << face.coord_mask<3>::operator std::string()<< std::endl;
                for(auto&& ind: face.indices) {
                    for(auto&& i: ind) {
                        std::cout << std::string(ccm.masked_vertex(i).mask()) << " ";
                    }
                }
                std::cout << std::endl;
            }
        }
    }


    return 0;
}


