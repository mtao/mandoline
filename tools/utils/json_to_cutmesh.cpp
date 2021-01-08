#include <fstream>
#include <filesystem>
#include <fmt/format.h>
// this set sould be removed after updating mtao core because prune vertices missed an include
#include <set>
#include <mtao/geometry/prune_vertices.hpp>
#include <mandoline/construction/remesh_self_intersections.hpp>
#include <mtao/json/bounding_box.hpp>
#include <mandoline/construction/construct.hpp>
#include "json_to_cutmesh.hpp"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
#include <igl/read_triangle_mesh.h>
#pragma GCC diagnostic pop

mandoline::CutCellMesh<3> json_to_cutmesh(const std::string& filename) {


    std::filesystem::path path(filename);


    std::ifstream ifs(filename);
    nlohmann::json js;
    ifs >> js;


    return json_to_cutmesh(js,path.parent_path());
}
mandoline::CutCellMesh<3> json_to_cutmesh(const nlohmann::json& js, const std::filesystem::path& relative_path) {
    auto bbox = js["grid"]["bounding_box"].get<Eigen::AlignedBox<double,3>>();
    std::array<int,3> shape;
    mtao::eigen::stl2eigen(shape) = mtao::json::json2vector<int,3>(js["grid"]["shape"]);






    auto mesh_path = relative_path / js["mesh"].get<std::string>();
    if(!std::filesystem::is_regular_file(mesh_path)) {
        // TODO: I should set some good enums for checking error codes someday right?
        //throw std::filesystem::filesystem_error("File not found", mesh_path,std::error_code(1,std::generic_category{}));
    }
    mtao::ColVecs3d V;
    mtao::ColVecs3i F;
    {

        bool do_remesh = js.contains("perform_remeshing") && js["perform_remeshing"].get<bool>();
        Eigen::MatrixXd VV;
        Eigen::MatrixXi FF;
        igl::read_triangle_mesh(mesh_path,VV,FF);
        V = VV.transpose();
        F = FF.transpose();
        if(do_remesh) {
            //auto t = mtao::logging::profiler("remesh_time",false,"remesh_profiler");
            std::tie(V,F) = mtao::geometry::prune(V,F,0);
            std::tie(V,F) = mandoline::construction::remesh_self_intersections(V,F);
        }
    }

    return mandoline::construction::from_bbox(V,F,bbox,shape);






}
