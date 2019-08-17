#include <mtao/geometry/grid/grid.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <mtao/geometry/grid/triangulation.hpp>
#include <igl/read_triangle_mesh.h>
#include <mtao/geometry/bounding_box.hpp>
#include <igl/writeOBJ.h>
#include <mtao/logging/profiler.hpp>



int main(int argc, char * argv[]) {
    if(argc < 3) {
        std::cout << "Usage: <scriptname> <input mesh> <grid edge>" << std::endl;
        return 1;
    }

    Eigen::MatrixXd iV,V;
    Eigen::MatrixXi iF,F;
    igl::read_triangle_mesh(argv[1], iV, iF);

    int N = std::atoi(argv[2]);
    std::array<int,3> NN{{N,N,N}};

    auto bbox = mtao::geometry::bounding_box(mtao::ColVecs3d(iV.transpose()));

    double bbox_offset = .1;

    mtao::Vec3d C = (bbox.min() + bbox.max()) / 2;
    mtao::Vec3d s = bbox.sizes() / 2;
    bbox.min() = C - (1 + bbox_offset) * s;
    bbox.max() = C + (1 + bbox_offset) * s;

    auto g = mtao::geometry::grid::Grid3d::from_bbox(bbox,NN,true);
    mtao::geometry::grid::GridTriangulator<std::decay_t<decltype(g)>> gt(g);
    Eigen::MatrixXd gV = gt.vertices().transpose();
    Eigen::MatrixXi gF = gt.faces().transpose();


    auto&& log_remesh = make_file_logger("MA_profiler","MA_profiler.log",mtao::logging::Level::Fatal,true);

    {
        auto t = mtao::logging::profiler("ma_time",false,"MA_profiler");

        igl::copyleft::cgal::mesh_boolean(iV,iF,gV,gF,igl::MESH_BOOLEAN_TYPE_RESOLVE,V,F);
    }

        auto&& dur = mtao::logging::profiler::durations();
        for(auto&& [pr,times]: dur) {
            auto&& [name,level] = pr;
            static const std::string pname = "MA_profiler"; if(name == pname) {
                static const std::string rmt = "ma_time";
                auto get_time = [&](const std::string& str) -> int {
                    auto dms = std::chrono::duration_cast<std::chrono::milliseconds>(times.at(str).first);
                    auto time = dms.count();
                    return time;
                };
                log_remesh.write(mtao::logging::Level::Fatal,  argv[1], ", ", N , ", " , get_time(rmt));
            }
        }

    igl::writeOBJ("output.obj",V,F);
    return 0;

}
