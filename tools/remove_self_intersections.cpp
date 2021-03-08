#include <igl/writeOBJ.h>

#include <iostream>
#include <mandoline/construction/tools/read_mesh.hpp>
#include <mandoline/construction/tools/remesh_self_intersections.hpp>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/geometry/prune_vertices.hpp>
#include <mtao/logging/profiler.hpp>
#include <set>

using namespace mandoline::construction;

int main(int argc, char* argv[]) {
    std::string filename = argv[1];

    auto [V, F] = tools::read_mesh(filename, false);
    // bool ret = igl::piecewise_constant_winding_number(F);
    auto&& log_remesh = make_file_logger("mesh_profiler", "mesh_profiler.log",
                                         mtao::logging::Level::Fatal, true);

    {
        auto t = mtao::logging::profiler("remesh_time", false, "mesh_profiler");
        std::tie(V, F) = mtao::geometry::prune(V, F, 0);
        std::tie(V, F) = tools::remesh_self_intersections(V, F);
    }
    auto&& dur = mtao::logging::profiler::durations();
    for (auto&& [pr, times] : dur) {
        auto&& [name, level] = pr;
        static const std::string pname = "mesh_profiler";
        if (name == pname) {
            static const std::string rmt = "remesh_time";
            auto get_time = [&](const std::string& str) -> int {
                auto dms =
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        times.at(str).first);
                auto time = dms.count();
                return time;
            };
            log_remesh.write(mtao::logging::Level::Fatal, argv[1], ", ",
                             get_time(rmt));
        }
    }
    std::set<std::array<int, 3>> Fs;
    for (int i = 0; i < F.cols(); ++i) {
        auto f = F.col(i);
        std::array<int, 3> v{{f(0), f(1), f(2)}};
        std::sort(v.begin(), v.end());
        Fs.emplace(v);
    }
    F = mtao::eigen::stl2eigen(Fs);

    igl::writeOBJ(argv[2], V.transpose(), F.transpose());

    return 0;
}
