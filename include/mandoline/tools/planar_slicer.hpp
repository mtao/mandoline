#pragma once

#include <balsa/eigen/types.hpp>
#include <Eigen/Geometry>
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mandoline/construction/cutdata.hpp>

namespace mandoline::tools {
class SliceGenerator : public mtao::geometry::grid::StaggeredGrid<double, 3> {
  public:
    using Base = mtao::geometry::grid::StaggeredGrid<double, 3>;
    SliceGenerator() {}
    SliceGenerator(const balsa::eigen::ColVecs3d &V, const balsa::eigen::ColVecs3i &F);
    std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i> slice(const balsa::eigen::Vec3d &origin, const balsa::eigen::Vec3d &direction);
    static Eigen::Affine3d get_transform(const balsa::eigen::Vec3d &origin, const balsa::eigen::Vec3d &direction);
    std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i> slice(const Eigen::Affine3d &t);
    void set_vertices(const balsa::eigen::ColVecs3d &V);
    Eigen::SparseMatrix<double> barycentric_map() const;

  private:
    void update_embedding(const balsa::eigen::ColVecs3d &V);
    balsa::eigen::ColVecs3d V;
    construction::CutData<3> data;
};

template<typename... Args>
std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i> slice(const balsa::eigen::ColVecs3d &V, const balsa::eigen::ColVecs3i &F, Args &&... args) {
    SliceGenerator sg(V, F);
    return sg.slice(std::forward<Args>(args)...);
}


}// namespace mandoline::tools
