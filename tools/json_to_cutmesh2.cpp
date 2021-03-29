#include "mandoline/construction/tools/json_to_cutmesh2.hpp"

#include <spdlog/spdlog.h>

#include <nlohmann/json.hpp>

int main(int argc, char* argv[]) {
    std::string filename = argv[1];

    std::filesystem::path path(filename);

    std::ifstream ifs(filename);
    nlohmann::json js;
    ifs >> js;

    const auto ppath = path.parent_path();

    auto ccm = mandoline::construction::tools::json_to_cutmesh2(js, ppath);

    std::ofstream ofs(argv[2]);

    auto V = ccm.vertices();
    auto E = ccm.edges();
    auto EE = E.array() + 1;
    for (int j = 0; j < V.cols(); ++j) {
        ofs << "v " << V.col(j).transpose() << " 0\n";
    }
    for (int j = 0; j < E.cols(); ++j) {
        ofs << "l " << EE.col(j).transpose() << "\n";
    }

    /*
    if (argc > 2) {
        ccm.write(argv[2]);
        return 0;
    } else if (js.contains("output_path")) {
        std::string opath = js["output_path"];
        if (!opath.empty()) {
            if (opath['.']) {
                ccm.write(ppath / opath);
                return 0;
            } else {
                ccm.write(opath);
                return 0;
            }
        } else {
            spdlog::error(
                "Must provide an output path, either through \"output_path\" "
                "in json or as a commandline argument");
            return 1;
        }
    }
    */

    return 0;
}
