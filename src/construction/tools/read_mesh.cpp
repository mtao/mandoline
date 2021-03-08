#include "mandoline/construction/tools/read_mesh.hpp"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
#include <igl/read_triangle_mesh.h>
#include <spdlog/spdlog.h>

#include <mtao/geometry/mesh/read_obj.hpp>
#pragma GCC diagnostic pop

#include "mandoline/construction/tools/preprocess_mesh.hpp"
namespace mandoline::construction::tools {

std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> read_mesh(
    const std::filesystem::path& path, bool preprocess) {
    mtao::ColVecs3d V;
    mtao::ColVecs3i F;
    if (path.filename().string().ends_with("obj")) {
        std::tie(V, F) = mtao::geometry::mesh::read_objD(path);
    } else {
        Eigen::MatrixXd VV;
        Eigen::MatrixXi FF;
        igl::read_triangle_mesh(path, VV, FF);
        V = VV.transpose();
        F = FF.transpose();
    }

    if (preprocess) {
        spdlog::info("Preprocessing mesh: {} verts, {} faces", V.cols(), F.cols());
        return preprocess_mesh(V, F);
    } else {
        spdlog::info("Returning mesh without preprocessing: {} verts, {} faces", V.cols(), F.cols());
        return {V, F};
    }
}
}  // namespace mandoline::construction::tools
