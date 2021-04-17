#include <cxxopts.hpp>
#include <mandoline/tools/plcurve_io.hpp>
#include <mandoline/construction/generator2.hpp>


int main(int argc, char *argv[]) {
    cxxopts::Options options("boundary_curve_to_cutmesh2", "A brief description");

    // clang-format off
    options.add_options()
        ("input", "input", cxxopts::value<std::string>())
        ("output", "output mesh file", cxxopts::value<std::string>())
        ("Nx", "Grid dimension in x", cxxopts::value<int>()->default_value("5"))
        ("Ny", "Grid dimension in y", cxxopts::value<int>()->default_value("5"))
        ("mx", "min value in x", cxxopts::value<double>()->default_value("0"))
        ("my", "min value in y", cxxopts::value<double>()->default_value("0"))
        ("Mx", "max value in x", cxxopts::value<double>()->default_value("10"))
        ("My", "max value in y", cxxopts::value<double>()->default_value("10"))
        ("h,help", "Print usage");
    // clang-format on
    options.parse_positional({ "input", "output" });
    options.positional_help({ "<input file> <output_file>" });

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }
    Eigen::AlignedBox<double, 2> bbox;
    bbox.min()
      << result["mx"].as<double>(),
      result["my"].as<double>();
    bbox.max()
      << result["Mx"].as<double>(),
      result["My"].as<double>();

    int Nx = result["Nx"].as<int>();
    int Ny = result["Ny"].as<int>();
    auto grid = mtao::geometry::grid::StaggeredGrid2d::from_bbox(bbox, { { Nx, Ny } });
    auto [oV, oE] = mandoline::tools::read_plcurve(result["input"].as<std::string>());


    mandoline::construction::CutCellGenerator<2> ccg(oV, grid);
    ccg.add_boundary_elements(oE.topRows<2>());
    ccg.bake();

    auto ccm = ccg.generate();
    mandoline::tools::write_cutmesh2_plcurve(ccm, result["output"].as<std::string>(), oE);
}

// OUTPUT:
// Vertices
// vertex_idx, xcoord, ycoord
//
// Edges
// edge_idx, vertex_idx1, vertex_idx2, left_element_idx, right_element_idx
//
// Elements
// element_idx, num_vertices, vertex_idx1, vertex_idx2, ..., vertex_idxn, componentID
