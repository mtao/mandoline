#pragma once
#include "mandoline/mesh3.hpp"
#include <balsa/eigen/types.hpp>


namespace mandoline::tools {
class MeshExploder {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    MeshExploder() {}
    MeshExploder(const CutCellMesh<3> &ccm);
    balsa::eigen::ColVecs3d V(double scale = 1.1, const std::set<int> &used_regions = {}) const;
    balsa::eigen::ColVecs3i F(const std::set<int> &used_regions = {}) const;
    balsa::eigen::ColVecs4d colors(const balsa::eigen::ColVecs4d &cell_colors, const std::set<int> &used_regions = {}) const;

    std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i> mesh(double scale = 1.1, const std::set<int> &used_regions = {}) const;
    std::tuple<balsa::eigen::ColVecs3d, balsa::eigen::ColVecs3i> mesh(size_t index, double scale = 1.1) const;
    size_t size() const { return Vs.size(); }
    balsa::eigen::ColVecs3d V(size_t index, double scale = 1.1) const;
    balsa::eigen::ColVecs3i F(size_t index) const;
    //region_scale in [0,1] interpolates the center between 0 (center is the grid center) and 1 (the center is the centroid of the region)
    void setCenters(double region_scale = 0.0);

    static bool valid_region(int region, const std::set<int> &used_regions = {}) {
        return used_regions.empty() || used_regions.find(region) != used_regions.end();
    }

  private:
    std::vector<balsa::eigen::ColVecs3d> Vs;
    std::vector<balsa::eigen::ColVecs3i> Fs;
    std::vector<int> regions;

    mtao::vector<balsa::eigen::Vec3d> RCs;
    mtao::vector<balsa::eigen::Vec3d> GCs;
    mtao::vector<balsa::eigen::Vec3d> Cs;
    balsa::eigen::Vec3d O;


    std::vector<int> offsets(const std::set<int> &used_regions = {}) const;
};
}// namespace mandoline::tools
