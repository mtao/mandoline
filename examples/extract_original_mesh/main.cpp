#include "mandoline/mesh3.hpp"
#include <mtao/cmdline_parser.hpp>
using namespace mtao::logging;



int main(int argc, char * argv[]) {

    mtao::CommandLineParser clp;
    clp.parse(argc, argv);

    if(clp.args().size() < 1) {
        fatal() << "No input mesh filename!";
        return {};
    } else if(clp.args().size() < 2) {
        fatal() << "No output obj filename!";
        return {};
    }

    std::string out_obj = clp.arg(1);


    std::ofstream ofs(out_obj);
    std::string input_cutmesh = clp.arg(0);
    auto ccm = mandoline::CutCellMesh<3>::from_file(input_cutmesh);
    auto& V = ccm.origV();
    auto& F = ccm.origF();
    for(int i = 0; i < V.cols(); ++i) {
        ofs << "v " << V.col(i).transpose() << std::endl;
    }
    for(int i = 0; i < F.cols(); ++i) {
        ofs << "f " << (F.col(i).array()+1).transpose() << std::endl;
    }





    return 0;
}


