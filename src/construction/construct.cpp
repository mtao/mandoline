#include "mandoline/construction/generator.hpp"
#include "mandoline/construction/construct.hpp"

namespace mandoline::construction {

    CutCellMesh<3> from_bbox(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const Eigen::AlignedBox<double,3>& bbox, const std::array<int,3>& cell_shape, int level, std::optional<double> threshold) {

        auto sg = CutCellMesh<3>::StaggeredGrid::from_bbox(bbox,cell_shape,false);
        return from_grid(V,F,sg,level,threshold);
    }
    CutCellMesh<3> from_grid(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const mtao::geometry::grid::StaggeredGrid3d& grid, int level, std::optional<double> threshold) {

        mtao::vector<mtao::Vec3d> stlp(V.cols());
        for(auto&& [i,v]: mtao::iterator::enumerate(stlp)) {
            v = V.col(i);
        }
        CutCellGenerator<3> ccg(stlp,grid, {});

        ccg.adaptive = true;
        ccg.adaptive_level = level;
        {
            auto t = mtao::logging::profiler("generator_bake",false,"profiler");
            ccg.add_boundary_elements(F);
            ccg.bake();
        }
        return ccg.generate();
    }
    CutCellMesh<3> from_grid(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const std::array<int,3>& cell_shape, int level, std::optional<double> threshold) {
        using Vec = mtao::Vec3d;
        auto sg = CutCellMesh<3>::StaggeredGrid(cell_shape, Vec::Ones());
        return from_grid(V,F,sg,level,threshold);
    }
}
