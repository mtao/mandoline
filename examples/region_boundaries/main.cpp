#include <cxxopts.hpp>
#include <mtao/geometry/mesh/write_obj.hpp>

#include "mandoline/operators/region_boundaries3.hpp"

int main(int argc, char* argv[]) {
    auto&& log = make_logger("profiler", mtao::logging::Level::All);

    cxxopts::Options clp(
        "region_boundaries",
        "Output the triangulated boundaries of a region of a mesh as an obj");

    // clang-format off
    clp.add_options()("r,region", "write per-region objs",
                      cxxopts::value<int>()->default_value("-1"))
        ("input_file", "input cutmesh file",cxxopts::value<std::string>())
        ("output_file", "output obj file",cxxopts::value<std::string>())
        ("h,help",
                                                                  "show help");

    // clang-format on
    clp.parse_positional({"input_file", "output_file"});

    auto result = clp.parse(argc, argv);

    int region = result["region"].as<int>();
    std::string input_file = result["input_file"].as<std::string>();
    std::string output_file = result["output_file"].as<std::string>();
    auto ccm = mandoline::CutCellMesh<3>::from_file(input_file);
    ccm.triangulate_faces();
    auto V = ccm.vertices();
    mtao::ColVecs3i F;
    std::vector<mtao::ColVecs3d> vertices;
    std::vector<mtao::ColVecs3i> faces;
    if (region < 0) {
        auto boundary_indices = mandoline::operators::region_boundaries(ccm);
        vertices.reserve(boundary_indices.size());
        faces.reserve(boundary_indices.size());
        for (auto&& index : boundary_indices) {
            std::tie(vertices.emplace_back(), faces.emplace_back()) =
                ccm.triangulate_face(index);
        }

    } else {
        auto boundary_map =
            mandoline::operators::region_boundaries(ccm, region);

        int vertex_offset = V.cols();

        for (auto&& [face, sgn] : boundary_map) {
            auto& v = vertices.emplace_back();
            auto& f = faces.emplace_back();
            std::tie(v, f) = ccm.triangulate_face(face);
            int vcols = v.cols();
            if (vcols > 0) {
                f.array() += vertex_offset;
                vertex_offset += vcols;
            }
            if (sgn) {
                mtao::RowVecXi tmp = f.row(1);
                f.row(1) = f.row(0);
                f.row(0) = tmp;
            }
        }
    }
    V = mtao::eigen::hstack(
        V, mtao::eigen::hstack_iter(vertices.begin(), vertices.end()));
    F = mtao::eigen::hstack_iter(faces.begin(), faces.end());

    mtao::geometry::mesh::write_objD(V, F, output_file);
    return 0;
}

