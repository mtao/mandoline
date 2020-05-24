#include "mandoline/construction/generator2.hpp"
#include "mandoline/construction/construct2.hpp"
#include <mtao/geometry/bounding_box.hpp>

namespace mandoline::construction {

CutCellMesh<2> from_bbox(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const Eigen::AlignedBox<double, 2> &bbox, const std::array<int, 2> &cell_shape, std::optional<double> threshold) {

    auto sg = CutCellMesh<2>::StaggeredGrid::from_bbox(bbox, cell_shape, false);
    return from_grid(V, F, sg, threshold);
}
CutCellMesh<2> from_grid(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const mtao::geometry::grid::StaggeredGrid2d &grid, std::optional<double> threshold) {

    mtao::vector<mtao::Vec2d> stlp(V.cols());
    for (auto &&[i, v] : mtao::iterator::enumerate(stlp)) {
        v = V.col(i);
    }
    CutCellGenerator<2> ccg(stlp, grid, {});

    {
        auto t = mtao::logging::profiler("generator_bake", false, "profiler");
        ccg.add_boundary_elements(F);
        ccg.bake();
    }
    return ccg.generate();
}
CutCellMesh<2> from_grid_unnormalized(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const std::array<int, 2> &cell_shape, std::optional<double> threshold) {
    using Vec = mtao::Vec2d;
    auto sg = CutCellMesh<2>::GridType(cell_shape, Vec::Ones());
    return from_grid(V, F, sg, threshold);
}

DeformingGeometryConstructor2::DeformingGeometryConstructor2() : _ccg{ std::make_unique<CutCellGenerator<2>>(std::array<int, 2>{ { 2, 2 } }) } {}

DeformingGeometryConstructor2::DeformingGeometryConstructor2(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const mtao::geometry::grid::StaggeredGrid2d &grid, std::optional<double> threshold) : _ccg(std::make_unique<CutCellGenerator<2>>(V, grid, threshold)) {

    assert(F.size() > 0);
    assert(V.size() > 0);
    {
        auto t = mtao::logging::profiler("generator_bake", false, "profiler");
        _ccg->add_boundary_elements(F);
        _ccg->bake();
        _dirty = false;
    }
}
DeformingGeometryConstructor2::~DeformingGeometryConstructor2() {}

void DeformingGeometryConstructor2::update_vertices(const mtao::ColVecs2d &V, const std::optional<double> &threshold) {
    _ccg->update_vertices(V, threshold);
    _dirty = true;
}
void DeformingGeometryConstructor2::update_topology(const mtao::ColVecs2i &F) {
    _ccg->set_boundary_elements(F);
}
void DeformingGeometryConstructor2::update_mesh(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const std::optional<double> &threshold) {

    update_vertices(V, threshold);
    update_topology(F);
    _dirty = true;
}
void DeformingGeometryConstructor2::update_grid(const mtao::geometry::grid::StaggeredGrid2d &g) {
    _ccg->update_grid(g);
    _dirty = true;
}
void DeformingGeometryConstructor2::set_vertices(const mtao::ColVecs2d &V, const std::optional<double> &threshold) {
    auto s = _ccg->vertex_shape();
    _ccg = std::make_unique<CutCellGenerator<2>>(V, *_ccg, threshold);
}
void DeformingGeometryConstructor2::set_mesh(const mtao::ColVecs2d &V, const mtao::ColVecs2i &E, const std::optional<double> &threshold) {
    set_vertices(V, threshold);
    update_topology(E);
    _dirty = true;
}
void DeformingGeometryConstructor2::bake() {
    if(!valid()) {
        spdlog::warn("DeformingGeometryConstructor2 can't bake when invalid");
    } else if (_dirty) {
        _ccg->clear();
        _ccg->bake();
        _dirty = false;
        spdlog::trace("DeoformingGeometryConstructor Done Baking");
    } else {
        spdlog::info("DeformingGeometryConstructor2: Redundant bake, nothing has changed");
    }
}
CutCellMesh<2> DeformingGeometryConstructor2::emit() const {
    return _ccg->generate();
}

DeformingGeometryConstructor2::DeformingGeometryConstructor2(DeformingGeometryConstructor2 &&o) : _ccg(std::move(o._ccg)), _dirty(o._dirty) {}
DeformingGeometryConstructor2 &DeformingGeometryConstructor2::operator=(DeformingGeometryConstructor2 &&o) {
    _ccg = std::move(o._ccg);
    _dirty = o._dirty;
    return *this;
}
bool DeformingGeometryConstructor2::valid() const {

    if (!_ccg) {
        spdlog::warn("No CCG object, constructor should have made it; ask a dev to fix this");
        return false;
    }
    if (!_ccg->valid()) {
        return false;
    }
    return true;
}
}// namespace mandoline::construction
