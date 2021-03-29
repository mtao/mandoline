#include "mandoline/construction/tools/json_to_cutmesh.hpp"

#include <spdlog/spdlog.h>

#include <nlohmann/json.hpp>

int main(int argc, char* argv[]) {
    std::string filename = argv[1];

    std::filesystem::path path(filename);

    std::ifstream ifs(filename);
    nlohmann::json js;
    ifs >> js;

    const auto ppath = path.parent_path();

    auto ccm = mandoline::construction::tools::json_to_cutmesh(js, ppath);
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

    return 0;
}
