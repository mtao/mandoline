#include <fmt/os.h>
#include <mandoline/operators/cell_indices.hpp>
#include <spdlog/spdlog.h>

#include <cxxopts.hpp>
#include <filesystem>
#include <fstream>
#include <mandoline/construction/tools/json_to_cutmesh.hpp>
#include <nlohmann/json.hpp>
using namespace mtao::logging;

using clock_type = std::chrono::high_resolution_clock;
int main(int argc, char* argv[]) {
    mtao::logging::make_logger("default", mtao::logging::Level::Off)
        .set_level(mtao::logging::Level::Off);
    mtao::logging::make_logger("profiler", mtao::logging::Level::Off)
        .set_level(mtao::logging::Level::Off);
    cxxopts::Options clp(
        "region_boundaries",
        "Output the triangulated boundaries of a region of a mesh as an obj");

    // clang-format off
    clp.add_options()
        ("input_json", "input json file",cxxopts::value<std::string>())
        ("output_csv", "output csv file",cxxopts::value<std::string>())
        ("s,samples", "samples per grid cell",cxxopts::value<int>()->default_value("10"))
        ("h,help","show help");

    // clang-format on

    clp.parse_positional({"input_json", "output_csv"});

    auto result = clp.parse(argc, argv);

    std::string input_file = result["input_json"].as<std::string>();
    std::string output_file = result["output_csv"].as<std::string>();
    int samples_per_cell = result["samples"].as<int>();

    std::filesystem::path path(input_file);
    if (!std::filesystem::exists(input_file)) {
        spdlog::error("mandoline: JSON file doesnt exist");
    }

    std::ifstream ifs(input_file);
    nlohmann::json js;
    ifs >> js;

    std::ofstream ofs(output_file);
    for (int N = 5; N < 500; N += 5) {
        js["grid"]["shape"] = {N, N, N};

        auto ccm = mandoline::construction::tools::json_to_cutmesh(
            js, path.parent_path());
        auto bb = ccm.bbox();
        mtao::ColVecs3d samples(3, N * N * N * samples_per_cell);
        samples.setRandom();
        samples = bb.sizes().eval().asDiagonal() * samples;
        samples.colwise() += bb.min();

        double average_grid_access = 0;
        double average_ccm_access = 0;

        const auto& grid = ccm.StaggeredGrid::vertex_grid();
        {
            mtao::VecXi I(samples.cols());
            auto start = clock_type::now();
            for (int j = 0; j < samples.cols(); ++j) {
                auto s = samples.col(j);
                auto coord = std::get<0>(grid.coord(s));
                I(j) = grid.index(coord);
            }
            auto end = clock_type::now();
            average_grid_access =
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      start)
                    .count() / double(samples.cols());
        }

        ccm.cache_vertices();
        {
            auto com = ccm.exterior_grid().cell_ownership_grid();
            mtao::VecXi I(samples.cols());
            auto start = clock_type::now();
            for (int j = 0; j < samples.cols(); ++j) {
                auto s = samples.col(j);
                I(j) = mandoline::operators::get_cell_index(ccm,com,s);
            }
            auto end = clock_type::now();
            average_ccm_access =
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      start)
                    .count() / double(samples.cols());
        }

        fmt::print("{},{},{},{},{}\n", N, ccm.num_cells(), ccm.num_cut_cells(),
                   average_grid_access, average_ccm_access);
        ofs << fmt::format("{},{},{},{},{}", N, ccm.num_cells(),
                           ccm.num_cut_cells(), average_grid_access,
                           average_ccm_access)
            << std::endl;
    }

    return 0;
}

