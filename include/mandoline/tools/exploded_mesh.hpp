#pragma once
#include "mandoline/mesh3.hpp"
#include <mtao/types.hpp>


namespace mandoline::tools {
    class MeshExploder {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
            MeshExploder() {}
            MeshExploder(const CutCellMesh<3>& ccm);
            mtao::ColVecs3d V(double scale=1.1, const std::set<int>& used_regions = {}) const;
            mtao::ColVecs3i F(const std::set<int>& used_regions = {}) const;

            std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> mesh(double scale=1.1, const std::set<int>& used_regions = {}) const;
            std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> mesh(size_t index, double scale=1.1) const;
            size_t size() const { return Vs.size(); }
            mtao::ColVecs3d V(size_t index, double scale=1.1) const;
            mtao::ColVecs3i F(size_t index) const;

            static bool valid_region(int region, const std::set<int>& used_regions = {}) {
                return used_regions.empty() ||  used_regions.find(region) != used_regions.end();
            }
        private:
            std::vector<mtao::ColVecs3d> Vs;
            std::vector<mtao::ColVecs3i> Fs;
            std::vector<int> regions;

            mtao::vector<mtao::Vec3d> Cs;
            mtao::Vec3d O;


            std::vector<int> offsets(const std::set<int>& used_regions = {}) const;
    };
}
