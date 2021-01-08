#include <Corrade/Utility/Arguments.h>
#include <fmt/format.h>

#include <iostream>
#include <nlohmann/json.hpp>

#include "validation/cutmesh_validation.hpp"

int main(int argc, char* argv[]) {
    Corrade::Utility::Arguments args;
    args.addArgument("filename")
        .setHelp("filename", ".cutmesh file to validate")
        .addBooleanOption('a', "all")
        .setHelp("all", "Run all tests")
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

    mandoline::CutCellMesh<3> ccm =
        mandoline::CutCellMesh<3>::from_proto(args.value("filename"));
    ccm.triangulate_faces(true);
    bool do_json = args.isSet("json");
    if(do_json) {
        js["filename"] = args.value("filename");
    }

    bool run_all = args.isSet("all");

    if (run_all || args.isSet("pcwn")) {
        auto [is_pcwn, is_not_pcwn] = pcwn_count(ccm);
        if (do_json) {
            js["pcwn"]["is"] = is_pcwn;
            js["pcwn"]["not"] = is_not_pcwn;
        } else {
            fmt::print("{} ", is_not_pcwn == 0);
        }
    }
    if (run_all || args.isSet("face-utilization")) {
        bool val = faces_fully_utilized(ccm);
        if (do_json) {
            js["face-utilization"] = val;
        } else {
            fmt::print("{} ", val);
        }
    }
    if (run_all || args.isSet("grid-utilization")) {
        auto [bad_volumes, unused] = grid_cells_fully_utilized_count(ccm);
        if (do_json) {
            js["grid-utilization"] = {{"bad_volumes", bad_volumes},
                                      {"unused", unused}};
        } else {
            fmt::print("{} ", bad_volumes == 0 && unused == 0);
        }
    }
    if (run_all || args.isSet("regions")) {
        auto [regions, input_regions] = region_counts(ccm);
        if (do_json) {
            js["regions"] = {{"regions", regions},
                             {"input_regions", input_regions}};

        } else {
            fmt::print("{} ", regions == input_regions);
        }
    }

    if (run_all || args.isSet("volume")) {
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
