#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <fstream>
// this set sould be removed after updating mtao core because prune vertices
// missed an include
#include <mandoline/construction/construct2.hpp>
#include <mtao/geometry/mesh/read_obj.hpp>
#include <mtao/json/bounding_box.hpp>
#include <set>

#include "json_to_cutmesh2.hpp"

mandoline::CutCellMesh<2> json_to_cutmesh2(const std::string& filename) {
    std::filesystem::path path(filename);
    if (!std::filesystem::exists(filename)) {
        spdlog::error("mandoline: JSON file doesnt exist");
    }

    std::ifstream ifs(filename);
    nlohmann::json js;
    ifs >> js;

    return json_to_cutmesh2(js, path.parent_path());
}
mandoline::CutCellMesh<2> json_to_cutmesh2(
    const nlohmann::json& js, const std::filesystem::path& relative_path) {
    auto bbox = js["grid"]["bounding_box"].get<Eigen::AlignedBox<double, 2>>();
    std::array<int, 2> shape;
    mtao::eigen::stl2eigen(shape) =
        mtao::json::json2vector<int, 2>(js["grid"]["shape"]);

    const std::string local_mesh_path = js["mesh"].get<std::string>();
    spdlog::info("Local mesh path: {}", local_mesh_path);
    if (local_mesh_path.empty()) {
        spdlog::error("No mesh path provided");
        return {};
    }
    std::filesystem::path mesh_path;
    if (local_mesh_path[0] == '.') {
        mesh_path = relative_path / local_mesh_path;
        spdlog::info("Found a dot");
    } else {
        mesh_path = local_mesh_path;
        spdlog::info("No dot, {}", local_mesh_path, std::string(mesh_path));
    }
    spdlog::info("Loading mesh path {}", std::string(mesh_path));
    if (!std::filesystem::is_regular_file(mesh_path)) {
        // TODO: I should set some good enums for checking error codes someday
        // right?
        // throw std::filesystem::filesystem_error("File not found",
        // mesh_path,std::error_code(1,std::generic_category{}));
    }
    auto [V, E] =
         mtao::geometry::mesh::read_obj2D(mesh_path);

    return mandoline::construction::from_bbox(V, E, bbox, shape);
}
