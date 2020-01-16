#include "mandoline/construction/construct_imgui.hpp"
#include <mtao/geometry/bounding_box.hpp>
namespace mandoline::construction {

    CutmeshGenerator_Imgui::CutmeshGenerator_Imgui(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const mtao::geometry::grid::StaggeredGrid3d& grid, int adaptive_level, std::optional<double> threshold)
        : DeformingGeometryConstructor(V,F,grid,adaptive_level,threshold)
          , bbox(grid.bbox())
          ,N(grid.vertex_shape())
          , use_cube(false)
          , adaptive_level(adaptive_level)
          , threshold(threshold)
    {
    }
    CutmeshGenerator_Imgui CutmeshGenerator_Imgui::create(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, double bbox_scale, const std::array<int,3>& N, int adaptive_level, std::optional<double> threshold)
    {
            auto bbox = mtao::geometry::bounding_box(V);
            bbox = mtao::geometry::expand_bbox(bbox,bbox_scale);
            auto s = bbox.sizes().eval();

            for(int i = 0; i < 3; ++i) {
                if(s(i)< 1e-5) {
                    bbox.min()(i) -= 1e-2;
                    bbox.max()(i) += 1e-2;
                }
            }
            auto g = mtao::geometry::grid::StaggeredGrid3d::from_bbox(bbox,N);
            return CutmeshGenerator_Imgui(V,F,g,adaptive_level,threshold);


    }

    void CutmeshGenerator_Imgui::set_bbox(const Eigen::AlignedBox<float,3>& bbox_, float scale) {
        bbox = mtao::geometry::expand_bbox(bbox_,scale);

        update_grid();
    }

    bool CutmeshGenerator_Imgui::gui() {
        bool dirty = false;

        if(ImGui::InputInt3("N", N.data()))  {
            dirty = true;
        }
        if(ImGui::InputFloat3("min", bbox.min().data()))  {
            bbox.min() = (bbox.min().array() < bbox.max().array()).select(bbox.min(),bbox.max());
            dirty = true;
        }
        if(ImGui::InputFloat3("max", bbox.max().data()))  {
            bbox.max() = (bbox.min().array() > bbox.max().array()).select(bbox.min(),bbox.max());
            dirty = true;
        }
        if(ImGui::Checkbox("Cubes", &use_cube)) {
            dirty = true;
        }
        if(ImGui::InputInt("Adaptive level", &adaptive_level))  {

            set_adaptivity(adaptive_level);
            dirty = true;
        }
        if(dirty) {
            update_grid();
        }
        return dirty;
    }

    void CutmeshGenerator_Imgui::update_grid() { 
        auto sg = staggered_grid();
        std::cout << "DX: ";
        std::cout << sg.dx().transpose() << std::endl;
        auto n = sg.vertex_shape();
        std::cout << "N: ";
        std::copy(n.begin(),n.end(),std::ostream_iterator<int>(std::cout," "));
        std::cout << std::endl;
        update_grid(sg);
        }

    mtao::geometry::grid::StaggeredGrid3d
        CutmeshGenerator_Imgui::staggered_grid() const { 
            return mtao::geometry::grid::StaggeredGrid3d::from_bbox
                (bbox.cast<double>(), N, use_cube);
        }

    mandoline::CutCellMesh<3> CutmeshGenerator_Imgui::generate() {
        bake();
        return emit();
    }
    void CutmeshGenerator_Imgui::update_vertices_and_bbox(const mtao::ColVecs3d& V, double bbox_scale, const std::optional<double>& threshold) {
        update_vertices(V, threshold);

        bbox = mtao::geometry::bounding_box(V).cast<float>();
        bbox = mtao::geometry::expand_bbox(bbox,float(bbox_scale));
        auto s = bbox.sizes().eval();

        for(int i = 0; i < 3; ++i) {
            if(s(i)< 1e-5) {
                bbox.min()(i) -= 1e-2;
                bbox.max()(i) += 1e-2;
            }
        }
        auto g = mtao::geometry::grid::StaggeredGrid3d::from_bbox(bbox.cast<double>(),N);
        update_grid();
    }
    void CutmeshGenerator_Imgui::update_mesh_and_bbox(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, double scale, const std::optional<double>& threshold) {

        update_vertices_and_bbox(V,scale,threshold);
        update_topology(F);
    }
}
