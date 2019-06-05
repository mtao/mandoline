#pragma once
#include "mandoline/mesh3.hpp"
#include <mtao/types.hpp>


namespace mandoline::tools {
    class MeshExploder {
        public:
            MeshExploder() {}
            MeshExploder(const CutCellMesh<3>& ccm);
            mtao::ColVecs3d V(double scale=1.1) const;
            mtao::ColVecs3i F() const;

            std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> mesh(double scale=1.1) const;
            std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> mesh(size_t index, double scale=1.1) const;
            size_t size() const { return Vs.size(); }
            mtao::ColVecs3d V(size_t index, double scale=1.1) const;
            mtao::ColVecs3i F(size_t index) const;

        private:
            std::vector<mtao::ColVecs3d> Vs;
            std::vector<mtao::ColVecs3i> Fs;

            mtao::vector<mtao::Vec3d> Cs;


            std::vector<int> offsets() const;
    };
}
