#include "make_cutmesh_from_cmdline.hpp"



int main(int argc, char * argv[]) {

    auto&& log = make_logger("profiler",mtao::logging::Level::All);
    auto clp = make_cutmesh_clparser();
    auto res = clp.parse(argc, argv);
    bool help_out = res.count("help");
    if (help_out) {
        std::cout << clp.help() << std::endl;
        return 0;
    }
    if(!bool(res.count("output"))) {
        mtao::logging::fatal() << "No output filename!";
        return {};
    }

    auto ccm = make_cutmesh(res);
    std::string output_prefix = res["output"].as<std::string>();;
    ccm.write(output_prefix);

    return 0;
}


