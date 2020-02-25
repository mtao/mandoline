#include "mandoline/tools/exploded_mesh.hpp"

namespace mandoline::tools {
MeshExploder::MeshExploder(const CutCellMesh<3> &ccm) {
    auto V = ccm.vertices();
    auto bb = ccm.bbox();
    O = (bb.min() + bb.max()) / 2;
    int off = 0;
    Vs.resize(ccm.cell_size());
    Fs.resize(ccm.cell_size());
    Cs.resize(ccm.cell_size());
    RCs.resize(ccm.cell_size());
    GCs.resize(ccm.cell_size());
    regions = ccm.regions();
    auto region_centroids = ccm.region_centroids();

    for (auto &&[i, c, v, f, rc, gc] : mtao::iterator::enumerate(ccm.cells(), Vs, Fs, RCs, GCs)) {
        std::tie(v, f) = c.get_mesh(V, ccm.faces());
        rc = region_centroids.col(c.region);
        gc = ccm.cell_grid().vertex(c.grid_cell);
    }

    auto &&AG = ccm.adaptive_grid();

    for (auto &&[idx, cell] : AG.cells()) {
        auto F = AG.triangulated(idx);
        std::tie(Vs[idx], Fs[idx]) = mtao::geometry::mesh::compactify(V, F);
        GCs[idx] = Vs[idx].rowwise().mean();
    }
    for (auto &&[c, i] : ccm.adaptive_grid_regions()) {
        RCs[c] = region_centroids.col(i);
    }
    setCenters();
}

void MeshExploder::setCenters(double region_scale) {

    for (auto &&[c, rc, gc] : mtao::iterator::zip(Cs, RCs, GCs)) {
        c = (1 - region_scale) * gc + region_scale * rc;
    }
}
mtao::ColVecs4d MeshExploder::colors(const mtao::ColVecs4d &cell_colors, const std::set<int> &used_regions) const {
    auto offs = offsets(used_regions);
    mtao::ColVecs4d C(4, offs.back());
    for (int i = 0; i < size(); ++i) {
        int off = offs[i];
        int size = offs[i + 1] - off;
        if (size > 0) {
            C.block(0, off, 4, size).colwise() = cell_colors.col(i);
        }
    }

    return C;
}

std::vector<int> MeshExploder::offsets(const std::set<int> &used_regions) const {
    std::vector<int> ret;
    ret.reserve(Vs.size() + 1);
    int off = 0;
    ret.push_back(off);
    for (auto &&[v, r] : mtao::iterator::zip(Vs, regions)) {
        if (valid_region(r, used_regions)) {
            off += v.cols();
        }
        ret.push_back(off);
    }
    return ret;
}

std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> MeshExploder::mesh(double scale, const std::set<int> &used_regions) const {
    return { V(scale, used_regions), F(used_regions) };
}
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> MeshExploder::mesh(size_t index, double scale) const {
    return { V(index, scale), F(index) };
}
mtao::ColVecs3d MeshExploder::V(double scale, const std::set<int> &used_regions) const {

    std::vector<mtao::ColVecs3d> V2;
    for (auto &&[i, r] : mtao::iterator::enumerate(regions)) {
        if (valid_region(r, used_regions)) {
            V2.emplace_back(V(i, scale));
        } else {
            V2.emplace_back();
        }
    }
    auto R = mtao::eigen::hstack_iter(V2.begin(), V2.end());
    return R;
}
mtao::ColVecs3i MeshExploder::F(const std::set<int> &used_regions) const {

    std::vector<mtao::ColVecs3i> F2;
    for (auto &&[F, o, r] : mtao::iterator::zip(Fs, offsets(used_regions), regions)) {
        if (valid_region(r, used_regions)) {
            F2.emplace_back(F.array() + o);
        } else {
            F2.emplace_back();
        }
    }
    return mtao::eigen::hstack_iter(F2.begin(), F2.end());
}
mtao::ColVecs3d MeshExploder::V(size_t index, double scale) const {
    auto R = Vs[index];
    R.colwise() -= O;
    R.colwise() += (scale - 1) * Cs[index];

    R.colwise() += O;

    return R;
}
mtao::ColVecs3i MeshExploder::F(size_t index) const {
    return Fs[index];
}
}// namespace mandoline::tools
