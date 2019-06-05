#pragma once
#include "mandoline/mesh3.hpp"
#include <mtao/types.hpp>


namespace mandoline::tools {
    class MeshExploder {
        public:
            MeshExploder(const CutCellMesh<3>& ccm);
            mtao::ColVecs3d V(double scale=1.1) const;
            mtao::ColVecs3i F() const;

            std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> mesh(double scale=1.1) const;

        private:
            std::vector<mtao::ColVecs3d> Vs;
            std::vector<mtao::ColVecs3i> Fs;

            mtao::vector<mtao::Vec3d> Cs;


            std::vector<int> offsets() const;
    };
}
