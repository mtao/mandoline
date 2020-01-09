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
    CutCellMesh<3> from_grid_unnormalized(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const std::array<int,3>& cell_shape, int level, std::optional<double> threshold) {
        using Vec = mtao::Vec3d;
        auto sg = CutCellMesh<3>::GridType(cell_shape, Vec::Ones());
        return from_grid(V,F,sg,level,threshold);
    }


    DeformingGeometryConstructor::DeformingGeometryConstructor(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const mtao::geometry::grid::StaggeredGrid3d& grid, int adaptive_level, std::optional<double> threshold): _ccg(new CutCellGenerator<3>(V,grid, threshold)) {

        _ccg->adaptive = true;
        _ccg->adaptive_level = adaptive_level;
        {
            auto t = mtao::logging::profiler("generator_bake",false,"profiler");
            _ccg->add_boundary_elements(F);
            _ccg->bake();
            _dirty = false;
        }
    }
    DeformingGeometryConstructor::~DeformingGeometryConstructor() {
        delete _ccg;
    }
    void DeformingGeometryConstructor::set_adaptivity(int res) {
        _ccg->adaptive_level = res;
        _dirty = true;
    }
    void DeformingGeometryConstructor::update_vertices(const mtao::ColVecs3d& V, const std::optional<double>& threshold) {
        _ccg->update_vertices(V, threshold);
        _dirty = true;
    }
    void DeformingGeometryConstructor::update_grid(const mtao::geometry::grid::StaggeredGrid3d& g) {
        _ccg->update_grid(g);
        _dirty = true;
    }
    void DeformingGeometryConstructor::bake() {
        if(_dirty) {
            _ccg->clear();
            _ccg->bake();
            _dirty = false;
        }
    }
    CutCellMesh<3> DeformingGeometryConstructor::emit() const {

        /*
            auto&& C = _ccg->crossings();
            for(auto&& c: C) {
                for(int d = 0; d < 3; ++d)
                {
                if(c.vertex().coord[d] == 0) {
                    if(c.vertex().quot(d) < 1e-5) {
                    std::cout << "Low quot " << std::string(c) << std::endl;
                    }
                } else if(c.vertex().coord[d] < 0) {
                    std::cout << "WTF <0 ?" << std::string(c) << std::endl;
                }
            }
            }
            */
        return _ccg->generate();
    }

    DeformingGeometryConstructor::DeformingGeometryConstructor(DeformingGeometryConstructor&& o):
                _ccg(o._ccg), _dirty(o._dirty)
            {o._ccg = nullptr;}
            DeformingGeometryConstructor& DeformingGeometryConstructor::operator=(DeformingGeometryConstructor&& o)
            {
                _ccg = o._ccg;
                o._ccg = nullptr;

                _dirty = o._dirty;
                return *this;
            }
}
