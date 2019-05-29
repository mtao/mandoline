#include "make_cutmesh_from_cmdline.hpp"



int main(int argc, char * argv[]) {

    auto&& log = make_logger("profiler",mtao::logging::Level::All);
    auto clp = make_cutmesh_clparser();
    if(clp.parse(argc, argv)) {

        std::string output_prefix = clp.arg(1);
        auto ccm = make_cutmesh(clp);
        ccm.write(output_prefix);
    }

    return 0;
}


