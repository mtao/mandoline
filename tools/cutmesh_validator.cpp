#include <Corrade/Utility/Arguments.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <iostream>
#include <mandoline/construction/construct.hpp>
#include <mandoline/construction/tools/read_mesh.hpp>
#include <mtao/geometry/bounding_box.hpp>
#include <mtao/geometry/mesh/read_obj.hpp>
#include <nlohmann/json.hpp>

#include "utils/json_to_cutmesh.hpp"
#include "validation/cutmesh_validation.hpp"

int main(int argc, char* argv[]) {
    Corrade::Utility::Arguments args;
    args.addArgument("filename")
        .setHelp("filename",
                 ".cutmesh file to validate or .obj to mesh or .json to mesh")
        .addBooleanOption('a', "all")
        .setHelp("all", "Run all tests")
        .addOption("bounding-box-scale", "0.3")
        .setHelp("bounding-box-scale",
                 "When loading an OBJ scale the bounding box")
        .addOption("grid-size", "5")
        .setHelp("grid-size", "Sidelength of the grid when loading an OBJ")
        .addBooleanOption('j', "json")
        .setHelp("json", "Output a json file with detailed statistics")
        .addBooleanOption('p', "pcwn")
        .setHelp("pcwn", "check cells for piecewise constant winding numbers")
        .addBooleanOption('f', "face-utilization")
        .setHelp("face-utilization",
                 "check for faces missing from the input mesh")
        .addBooleanOption('g', "grid-utilization")
        .setHelp("grid-utilization", "check for grid cells missing")
        .addBooleanOption('r', "regions")
        .setHelp("regions",
                 "check that the result has the same number as regions as the "
                 "input mesh")
        .addBooleanOption('v', "volume")
        .setHelp("volume", "Check if volume is equal to the input")
        .addBooleanOption('b', "boundary")
        .setHelp("boundary",
                 "Confirm that each face is the boundary of two cells or is a "
                 "grid boundary")
        .addBooleanOption('e', "external-valence")
        .setHelp("external-valence",
                 "Check that purely cubical cells have valence 6");

    args.parse(argc, argv);
    using json = nlohmann::json;
    json js;

    mandoline::CutCellMesh<3> ccm;
    const std::string filename = args.value("filename");
    if (filename.ends_with("obj")) {
        auto [V,F] = mandoline::construction::tools::read_mesh(filename,true);
        auto bb = mtao::geometry::bounding_box(V);
        mtao::Vec3d s = bb.sizes();
        s *= std::stof(args.value("bounding-box-scale"));
        spdlog::info("Bounding box: = {} => {}", fmt::join(bb.min(), ","),
                     fmt::join(bb.max(), ","));

        spdlog::info("s = {}", fmt::join(s, ","));
        bb.min() -= s;
        bb.max() += s;
        spdlog::info("Bounding box: = {} => {}", fmt::join(bb.min(), ","),
                     fmt::join(bb.max(), ","));
        int N = std::stoi(args.value("grid-size"));
        spdlog::info("N = {}", N);
        ccm = mandoline::construction::from_bbox(V, F, bb,
                                                 std::array<int, 3>{{N, N, N}});

    } else if (filename.ends_with("json")) {
        ccm = json_to_cutmesh(filename);
    } else {
        ccm = mandoline::CutCellMesh<3>::from_proto(filename);
    }
    auto tmp_path = std::filesystem::temp_directory_path();
    auto output_cutmesh = tmp_path / "validation_output.cutmesh";
    auto output_config = tmp_path / "cutmesh_config.json";
    spdlog::info(
            "Outputting generated mesh to /tmp/validation_output.cutmesh");
    ccm.write(output_cutmesh);
    ccm.triangulate_faces(true);
    bool do_json = args.isSet("json");
    if (do_json) {
        js["filename"] = args.value("filename");
    }

    bool run_all = args.isSet("all");

    if (run_all || args.isSet("pcwn")) {
        spdlog::info("Checking cells for PCWN");
        auto [is_pcwn, is_not_pcwn] = pcwn_count(ccm);
        if (do_json) {
            js["pcwn"]["is"] = is_pcwn;
            js["pcwn"]["not"] = is_not_pcwn;
        } else {
            fmt::print("{} ", is_not_pcwn == 0);
        }
    }
    if (run_all || args.isSet("face-utilization")) {
        spdlog::info("Checking face utilization");
        bool val = faces_fully_utilized(ccm);
        if (do_json) {
            js["face-utilization"] = val;
        } else {
            fmt::print("{} ", val);
        }
    }
    if (run_all || args.isSet("grid-utilization")) {
        spdlog::info("Checking grid utilization");
        auto [bad_volumes, unused] = grid_cells_fully_utilized_count(ccm);
        if (do_json) {
            js["grid-utilization"] = {{"bad_volumes", bad_volumes},
                                      {"unused", unused}};
        } else {
            fmt::print("{} ", bad_volumes == 0 && unused == 0);
        }
    }
    if (run_all || args.isSet("regions")) {
        spdlog::info("Checking regions");
        auto [regions, input_regions] = region_counts(ccm);
        if (do_json) {
            js["regions"] = {{"regions", regions},
                             {"input_regions", input_regions}};

        } else {
            fmt::print("{} ", regions == input_regions);
        }
    }

    if (run_all || args.isSet("volume")) {
        spdlog::info("Checking volume");
        bool val = volume_check(ccm);
        if (do_json) {
            js["volume"] = val;
        } else {
            fmt::print("{} ", val);
        }
    }

    if (run_all || args.isSet("boundary")) {
        bool val = paired_boundary(ccm);
        if (do_json) {
            fmt::print("{} ", val);

        } else {
            js["boundary"] = val;
        }
    }

    if (run_all || args.isSet("external-valence")) {
        bool val = exterior_cell_valence_counts(ccm);
        if (do_json) {
            js["external-valence"] = val;
        } else {
            fmt::print("{}", val);
        }
    }
    if (do_json) {
        std::cout << js.dump(2) << std::endl;
    } else {
        fmt::print("\n");
    }
}
