#include <cxxopts.hpp>
#include <mandoline/mesh3.hpp>


int main(int argc, char *argv[]) {
    cxxopts::Options options("add_triangulation_to_cutmesh", "Add a triangulation to a cutmesh");

    // clang-format off
    options.add_options()
        ("input", "input cutmesh", cxxopts::value<std::string>())
        ("output", "output cutmesh", cxxopts::value<std::string>())
        ("n,new-vertices", "Potentially add new vertices", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage");
    // clang-format on
    options.parse_positional({ "input", "output" });
    options.positional_help({ "<input file> <output_file>" });

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    std::string input = result["input"].as<std::string>();
    std::string output = result["output"].as<std::string>();
    bool new_vertices = result["new-vertices"].as<bool>();
    auto ccm = mandoline::CutCellMesh<3>::from_proto(input);

    ccm.triangulate_faces(new_vertices);

    ccm.write(output);

    return 0;

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
