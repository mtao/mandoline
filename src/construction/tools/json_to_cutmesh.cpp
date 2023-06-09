#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <fstream>
// this set sould be removed after updating mtao core because prune vertices
// missed an include
#include <mtao/json/bounding_box.hpp>
#include <set>

#include "mandoline/construction/construct.hpp"
#include "mandoline/construction/tools/json_to_cutmesh.hpp"
#include "mandoline/construction/tools/read_mesh.hpp"

namespace mandoline::construction::tools {

CutCellMesh<3> json_to_cutmesh(const std::string& filename) {
    std::filesystem::path path(filename);
    if (!std::filesystem::exists(filename)) {
        spdlog::error("mandoline: JSON file doesnt exist");
    }

    std::ifstream ifs(filename);
    nlohmann::json js;
    ifs >> js;

    return json_to_cutmesh(js, path.parent_path());
}
CutCellMesh<3> json_to_cutmesh(const nlohmann::json& js,
                               const std::filesystem::path& relative_path) {
    auto bbox = js["grid"]["bounding_box"].get<Eigen::AlignedBox<double, 3>>();
    std::array<int, 3> shape;
    mtao::eigen::stl2eigen(shape) =
        mtao::json::json2vector<int, 3>(js["grid"]["shape"]);

    const std::string local_mesh_path = js["mesh"].get<std::string>();
    spdlog::trace("Local mesh path: {}", local_mesh_path);
    if (local_mesh_path.empty()) {
        spdlog::error("No mesh path provided");
        return {};
    }
    std::filesystem::path mesh_path;
    if (local_mesh_path[0] == '.') {
        mesh_path = relative_path / local_mesh_path;
        spdlog::trace("Found a dot");
    } else {
        mesh_path = local_mesh_path;
        spdlog::trace("No dot, {}", local_mesh_path, std::string(mesh_path));
    }
    spdlog::trace("Loading mesh path {}", std::string(mesh_path));
    if (!std::filesystem::is_regular_file(mesh_path)) {
        // TODO: I should set some good enums for checking error codes someday
        // right?
        // throw std::filesystem::filesystem_error("File not found",
        // mesh_path,std::error_code(1,std::generic_category{}));
    }
    bool do_remesh =
        js.contains("perform_remeshing") && js["perform_remeshing"].get<bool>();
    auto [V, F] = read_mesh(mesh_path, do_remesh);

    spdlog::info("Grid shape [{}] {} => {}", fmt::join(shape, ","),
                 fmt::join(bbox.min(), ","), fmt::join(bbox.max(), ","));
    return from_bbox(V, F, bbox, shape);
}
}  // namespace mandoline::construction::tools