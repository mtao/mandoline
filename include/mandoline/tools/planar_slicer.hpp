#pragma once

#include <mtao/types.hpp>
#include <Eigen/Geometry>

namespace mandoline::tools {
    class SliceGenerator {
        public:
            SliceGenerator() {}
            SliceGenerator(const mtao::Vec3d& origin, const mtao::Vec3d& direction);
            SliceGenerator(const Eigen::Affine3d& t): transform(t) {}
            std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> slice(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F) const;
            //std::tuple<mtao::ColVecs3d& V, mtao::ColVecs3i& F> slice(const mtao::ColVecs3d& V, const CutData& cut_data) const;
        private:
            Eigen::Affine3d transform;
    };

    template <typename... Args>
    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> slice(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, Args&&... args) {
        SliceGenerator sg(std::forward<Args>(args)...);
        return sg.slice(V,F);
    }


}

