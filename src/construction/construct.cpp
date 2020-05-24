#include "mandoline/construction/generator3.hpp"
#include "mandoline/construction/construct.hpp"
#include <mtao/geometry/bounding_box.hpp>

namespace mandoline::construction {

CutCellMesh<3> from_bbox(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const Eigen::AlignedBox<double, 3> &bbox, const std::array<int, 3> &cell_shape, int level, std::optional<double> threshold) {

    auto sg = CutCellMesh<3>::StaggeredGrid::from_bbox(bbox, cell_shape, false);
    return from_grid(V, F, sg, level, threshold);
}
CutCellMesh<3> from_grid(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const mtao::geometry::grid::StaggeredGrid3d &grid, int level, std::optional<double> threshold) {

    mtao::vector<mtao::Vec3d> stlp(V.cols());
    for (auto &&[i, v] : mtao::iterator::enumerate(stlp)) {
        v = V.col(i);
    }
    CutCellGenerator<3> ccg(stlp, grid, {});

    ccg.adaptive_level = level;
    {
        auto t = mtao::logging::profiler("generator_bake", false, "profiler");
        ccg.add_boundary_elements(F);
        ccg.bake();
    }
    return ccg.generate();
}
CutCellMesh<3> from_grid_unnormalized(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const std::array<int, 3> &cell_shape, int level, std::optional<double> threshold) {
    using Vec = mtao::Vec3d;
    auto sg = CutCellMesh<3>::GridType(cell_shape, Vec::Ones());
    return from_grid(V, F, sg, level, threshold);
}


DeformingGeometryConstructor::DeformingGeometryConstructor(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const mtao::geometry::grid::StaggeredGrid3d &grid, int adaptive_level, std::optional<double> threshold) : _ccg(std::make_unique<CutCellGenerator<3>>(V, grid, threshold)) {

    _ccg->adaptive_level = adaptive_level;
    assert(F.size() > 0);
    assert(V.size() > 0);
    {
        auto t = mtao::logging::profiler("generator_bake", false, "profiler");
        _ccg->add_boundary_elements(F);
        _ccg->bake();
        _dirty = false;
    }
}
DeformingGeometryConstructor::DeformingGeometryConstructor() : _ccg{ std::make_unique<CutCellGenerator<3>>(std::array<int, 3>{ { 2, 2, 2 } }) } {}
DeformingGeometryConstructor::~DeformingGeometryConstructor() {}
void DeformingGeometryConstructor::set_adaptivity(int res) {
    _ccg->adaptive_level = res;
    _dirty = true;
}
void DeformingGeometryConstructor::update_vertices(const mtao::ColVecs3d &V, const std::optional<double> &threshold) {
    _ccg->update_vertices(V, threshold);
    _dirty = true;
}
void DeformingGeometryConstructor::update_topology(const mtao::ColVecs3i &F) {
    _ccg->set_boundary_elements(F);
}
void DeformingGeometryConstructor::set_vertices(const mtao::ColVecs3d &V, const std::optional<double> &threshold) {
    auto s = _ccg->vertex_shape();
    _ccg = std::make_unique<CutCellGenerator<3>>(V, *_ccg, threshold);
}
void DeformingGeometryConstructor::set_mesh(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const std::optional<double> &threshold) {
    set_vertices(V, threshold);
    update_topology(F);
    _dirty = true;
}
void DeformingGeometryConstructor::update_mesh(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const std::optional<double> &threshold) {

    update_vertices(V, threshold);
    update_topology(F);
    _dirty = true;
}
void DeformingGeometryConstructor::update_grid(const mtao::geometry::grid::StaggeredGrid3d &g) {
    _ccg->update_grid(g);
    _dirty = true;
}
void DeformingGeometryConstructor::bake() {
    if (_dirty && valid()) {
        _ccg->clear();
        _ccg->bake();
        _dirty = false;
    }
    spdlog::trace("DeoformingGeometryConstructor Done Baking");
}
CutCellMesh<3> DeformingGeometryConstructor::emit() const {
    return _ccg->generate();
}

DeformingGeometryConstructor::DeformingGeometryConstructor(DeformingGeometryConstructor &&o) : _ccg(std::move(o._ccg)), _dirty(o._dirty) {}
DeformingGeometryConstructor &DeformingGeometryConstructor::operator=(DeformingGeometryConstructor &&o) {
    _ccg = std::move(o._ccg);
    _dirty = o._dirty;
    return *this;
}
bool DeformingGeometryConstructor::valid() const {

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
