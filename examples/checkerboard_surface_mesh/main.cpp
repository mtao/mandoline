#include <fmt/core.h>
#include <spdlog/spdlog.h>
#include <tbb/parallel_for_each.h>

#include <cxxopts.hpp>
#include <mandoline/mesh3.hpp>
#include <mtao/filesystem/prepend_to_filename.hpp>
#include <mtao/geometry/mesh/write_obj.hpp>
#include <mtao/geometry/mesh/write_ply.hpp>
#include <mutex>
#include <vector>

using namespace mandoline;
using CoordType = std::array<int, 3>;
auto possible_cells(const CutCellMesh<3>& ccm, const std::vector<int>& face)
    -> std::set<CoordType> {
    if (face.empty()) {
        return {};
    }
    auto possible = [&](int idx) -> std::set<std::array<int, 3>> {
        if (ccm.is_grid_vertex(idx)) {
            return Vertex<3>(ccm.vertex_grid().unindex(idx)).possible_cells();
        } else {
            return {};
        }
    };

    std::set<CoordType> possibles = possible(face[0]);

    for (auto&& f : face) {
        auto s = possible(f);
        std::set<CoordType> i;
        std::set_intersection(possibles.begin(), possibles.end(), s.begin(),
                              s.end(), std::inserter(i, i.end()));
        possibles = std::move(i);
        if (possibles.empty()) {
            return {};
        }
    }

    return possibles;
}

struct Polys {
    // std::vector<mtao::ColVecs3d> vertices;
    std::vector<mtao::ColVecs3i> faces;
    Polys(const CutCellMesh<3>& ccm, auto&& filter) {
        // auto V = ccm.vertices();

        faces.reserve(ccm.mesh_cut_faces().size());
        // vertices.reserve(ccm.mesh_cut_faces().size());

        std::mutex mutex;
        int voff = 0;

        // for (auto&& [idx, _] : ccm.mesh_cut_faces()) {
        const auto& mcf = ccm.mesh_cut_faces();
        tbb::parallel_for_each(mcf.begin(), mcf.end(), [&](auto&& pr) {
            int idx = std::get<0>(pr);
            const auto& face = ccm.faces().at(idx);
            auto coord = face.get_min_coord(
                [&](int idx) { return ccm.masked_vertex(idx).coord; });
            bool use = filter(coord);

            if (use) {
                // auto [V, F] = ccm.compact_triangulated_face(idx);
                auto F = face.triangulate_fan();

                std::scoped_lock lock(mutex);
                // F.array() += voff;
                // voff += V.cols();

                // vertices.emplace_back(std::move(V));
                faces.emplace_back(std::move(F));
            }
        });
    }
    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> mesh(
        const CutCellMesh<3>& ccm) const {
        auto V = ccm.vertices();
        // auto V = mtao::eigen::hstack_iter(vertices.begin(), vertices.end());
        auto F = mtao::eigen::hstack_iter(faces.begin(), faces.end());
        return mtao::geometry::mesh::compactify(V, F);
    }
};

int main(int argc, char* argv[]) {
    auto&& log = make_logger("profiler", mtao::logging::Level::All);

    cxxopts::Options clp(
        "checkerboard_surface_mesh",
        "Shows the original mesh with a checkerboard coloring on teh surface");

    // clang-format off
    clp.add_options()("r,region", "write per-region objs",
                      cxxopts::value<int>()->default_value("-1"))
        ("input_file", "input cutmesh file",cxxopts::value<std::string>())
        ("output_file", "output ply file",cxxopts::value<std::string>())
        ("light_color", "lighter color of checkerboard",cxxopts::value<std::string>()->default_value("c18741"))
        ("dark_color", "darkercolor of checkerboard",cxxopts::value<std::string>()->default_value("9b6241"))
        ("o,obj","output two obj files, output is treated as a fmt format");
        ("h,help","show help");

    // clang-format on
    clp.parse_positional({"input_file", "output_file"});

    auto result = clp.parse(argc, argv);

    auto valid_hex_char = [](char c) -> bool {
        if ('0' <= c && c <= '9') {
            return true;
        } else if ('a' <= c && c <= 'f') {
            return true;
        } else if ('A' <= c && c <= 'F') {
            return true;
        } else {
            return false;
        }
    };
    auto hex_number = [&](char c) -> int {
        if ('0' <= c && c <= '9') {
            return c - '0';
        } else if ('a' <= c && c <= 'f') {
            return 10 + (c - 'a');
        } else if ('A' <= c && c <= 'F') {
            return 10 + (c - 'A');
        } else {
            return {};
        }
    };

    auto valid_hex = [&](const std::string& str) -> bool {
        if (str.size() != 6 && str.size() != 8) {
            spdlog::error("hex color [{}] needs 6 or 8 characters, got {}", str,
                          str.size());
            return false;
        }
        for (auto&& c : str) {
            if (!valid_hex_char(c)) {
                spdlog::error("hex color [{}] contained an invalid char [{}]",
                              str, c);
                return false;
            }
        }
        return true;
    };

    auto parse_hex = [&](const std::string& str) -> Eigen::Vector4d {
        Eigen::Vector4d r;
        if (!valid_hex(str)) {
            return r;
        }
        r.x() = (16 * hex_number(str[0]) + hex_number(str[1])) / 256.;
        r.y() = (16 * hex_number(str[2]) + hex_number(str[3])) / 256.;
        r.z() = (16 * hex_number(str[4]) + hex_number(str[5])) / 256.;
        if (str.size() == 8) {
            r.w() = (16 * hex_number(str[6]) + hex_number(str[7])) / 256.;
        } else {
            r.w() = 1.0;
        }
        return r;
    };

    const std::string light = result["light_color"].as<std::string>();
    const std::string dark = result["dark_color"].as<std::string>();

    const mtao::Vec4d light_col = parse_hex(light);
    const mtao::Vec4d dark_col = parse_hex(dark);

    const std::string input_file = result["input_file"].as<std::string>();
    const std::string output_file = result["output_file"].as<std::string>();

    const bool output_obj = result["obj"].count() > 0;

    auto ccm = mandoline::CutCellMesh<3>::from_file(input_file);
    mtao::ColVecs4d C;
    mtao::ColVecs3i F;
    std::vector<mtao::ColVecs3d> vertices;
    std::vector<mtao::ColVecs3i> faces;
    std::vector<mtao::ColVecs4d> colors;

    faces.reserve(ccm.mesh_cut_faces().size());
    // ccm.triangulate_faces(true);

    Polys lP(ccm, [](const CoordType& c) -> bool {
        int sum = std::accumulate(c.begin(), c.end(), int(0));
        return sum % 2 == 0;
    });
    Polys dP(ccm, [](const CoordType& c) -> bool {
        int sum = std::accumulate(c.begin(), c.end(), int(0));
        return sum % 2 != 0;
    });

    // V = mtao::eigen::hstack_iter(vertices.begin(), vertices.end());
    // C = mtao::eigen::hstack_iter(colors.begin(), colors.end());
    // F = mtao::eigen::hstack_iter(faces.begin(), faces.end());

    {
        auto [lV, lF] = lP.mesh(ccm);
        auto [dV, dF] = dP.mesh(ccm);

        if (output_obj) {
            std::string light_filename =
                fmt::vformat(output_file, fmt::make_format_args("light"));
            std::string dark_filename =
                fmt::vformat(output_file, fmt::make_format_args("dark"));
            if (light_filename == dark_filename) {
                spdlog::warn(
                    "Format was invalid [{}], both outputs had the same name "
                    "[{}]",
                    output_file, light_filename, dark_filename);
                light_filename =
                    mtao::filesystem::prepend_to_filename(input_file, "light_")
                        .replace_extension(".obj");
                dark_filename =
                    mtao::filesystem::prepend_to_filename(input_file, "dark_")
                        .replace_extension(".obj");
                spdlog::warn("Output  filenames set to {} and {}",
                             light_filename, dark_filename);
            }
            mtao::geometry::mesh::write_objD(lV, lF, light_filename);
            mtao::geometry::mesh::write_objD(dV, dF, dark_filename);
        } else {
            auto V = mtao::eigen::hstack(lV, dV);
            auto F = mtao::eigen::hstack(lF, dF.array() + lV.cols());
            auto C =
                mtao::eigen::hstack(light_col.rowwise().replicate(lV.cols()),
                                    dark_col.rowwise().replicate(dV.cols()));

            mtao::geometry::mesh::write_plyD(V, C.topRows<3>(), F, output_file);
        }
    }
    return 0;
}

