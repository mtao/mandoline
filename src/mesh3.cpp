#include "mandoline/mesh3.hpp"

#include <igl/winding_number.h>
#include <mtao/geometry/volume.h>
#include <spdlog/spdlog.h>

#include <mtao/geometry/bounding_box.hpp>
#include <mtao/geometry/kdtree.hpp>
#include <mtao/geometry/mesh/halfedge.hpp>
#include <mtao/geometry/mesh/separate_elements.hpp>
#include <mtao/geometry/prune_vertices.hpp>
#include <mtao/iterator/enumerate.hpp>
#include <mtao/iterator/range.hpp>
#include <mtao/logging/logger.hpp>

#include "mandoline/diffgeo_utils.hpp"
#include "mandoline/operators/boundary3.hpp"
#include "mandoline/operators/centroids3.hpp"
#include "mandoline/operators/interpolation3.hpp"
#include "mandoline/operators/masks.hpp"
#include "mandoline/operators/volume3.hpp"
#include "mandoline/proto_util.hpp"
#if defined(MTAO_TBB_ENABLED)
#include <tbb/parallel_for_each.h>
#endif

namespace mandoline {

size_t CutCellMesh<3>::cell_size() const {
    return m_exterior_grid.num_cells() + m_cells.size();
}
size_t CutCellMesh<3>::cut_face_size() const { return m_faces.size(); }
size_t CutCellMesh<3>::num_cells() const {
    return m_exterior_grid.num_cells() + num_cut_cells();
}
size_t CutCellMesh<3>::cut_cell_size() const { return m_cells.size(); }
size_t CutCellMesh<3>::num_cut_cells() const { return m_cells.size(); }
size_t CutCellMesh<3>::face_size() const {
    return m_exterior_grid.num_faces() + cut_face_size();
}

size_t CutCellMesh<3>::num_faces() const {
    return m_exterior_grid.num_faces() + cut_face_size();
}
size_t CutCellMesh<3>::num_cut_faces() const { return m_faces.size(); }

bool CutCellMesh<3>::is_cut_face(size_t index) const {
    return index < num_cut_faces();
}
bool CutCellMesh<3>::is_cut_cell(int index) const {
    return index >= 0 && index < m_cells.size();
}
bool CutCellMesh<3>::is_exterior_cell(int index) const {
    return index >= 0 && index >= m_cells.size();
}
bool CutCellMesh<3>::is_mesh_face(int idx) const {
    return m_mesh_cut_faces.find(idx) != m_mesh_cut_faces.end();
}
auto CutCellMesh<3>::cell_volumes() const -> VecX {
    return operators::cell_volumes(*this);
}
auto CutCellMesh<3>::face_centroids() const -> mtao::ColVecs3d {
    return operators::face_centroids(*this);
}
auto CutCellMesh<3>::cell_centroids() const -> mtao::ColVecs3d {
    return operators::cell_centroids(*this);
}
auto CutCellMesh<3>::dual_vertices() const -> ColVecs {
    return cell_centroids();
}
int CutCellMesh<3>::local_grid_cell_index(const VecCRef &p) const {
    auto [c, q] = StaggeredGrid::coord(p);
    if (!StaggeredGrid::cell_grid().valid_index(c)) {
        return -1;
    }

    return -1;
}
int CutCellMesh<3>::world_grid_cell_index(const VecCRef &p) const {
    return local_grid_cell_index(vertex_grid().local_coord(p));
}
std::vector<std::set<int>> CutCellMesh<3>::cell_faces() const {
    std::vector<std::set<int>> ret(m_cells.size());
    std::transform(m_cells.begin(), m_cells.end(), ret.begin(), [](auto &&V) {
        std::set<int> ret;
        std::transform(V.begin(), V.end(), std::inserter(ret, ret.end()),
                       [](auto &&v) { return std::get<0>(v); });

        return ret;
    });
    return ret;
}
std::set<int> CutCellMesh<3>::cell_faces(int index) const {
    auto &V = m_cells.at(index);
    std::set<int> ret;
    std::transform(V.begin(), V.end(), std::inserter(ret, ret.end()),
                   [](auto &&v) { return std::get<0>(v); });

    return ret;
}

void CutCellMesh<3>::write(const std::string &filename) const {
    GOOGLE_PROTOBUF_VERIFY_VERSION;
    std::stringstream ss;
    std::ofstream ofs(filename, std::ios::binary);
    protobuf::CutMeshProto cmp;
    serialize(cmp);
    cmp.SerializeToOstream(&ofs);
}

std::vector<bool> CutCellMesh<3>::boundary_faces() const {
    std::vector<bool> ret(m_faces.size(), false);
    for (auto &&[m, _] : m_mesh_cut_faces) {
        ret[m] = true;
    }
    return ret;
}

std::vector<int> CutCellMesh<3>::regions(bool boundary_sign_regions) const {
    std::vector<int> ret(cell_size(), 1);

    if (boundary_sign_regions) {
        for (auto &&[idx, c] : mtao::iterator::enumerate(m_cells)) {
            Edge counts{{0, 0}};  // 1 -1
            for (auto &&[f, b] : c) {
                if (is_mesh_face(f)) {
                    counts[b ? 0 : 1]++;
                }
            }
            if (counts[0] + counts[1] > 0) {
                if (counts[0] > counts[1]) {
                    ret[idx] = 2;
                } else {
                    ret[idx] = 3;
                }
            }
        }
    } else {
        for (int i = 0; i < m_cells.size(); ++i) {
            ret[i] = m_cells[i].region;
        }
        for (auto &&[c, r] : m_adaptive_grid_regions) {
            ret[c] = r;
        }
    }
    return ret;
}

std::vector<std::array<std::set<int>, 2>> CutCellMesh<3>::face_regions() const {
    auto R = regions();
    std::copy(R.begin(), R.end(), std::ostream_iterator<int>(std::cout, ","));
    std::vector<std::array<std::set<int>, 2>> ret(
        *std::max_element(R.begin(), R.end()) + 1);
    for (auto &&[cidx, c] : mtao::iterator::enumerate(cells())) {
        for (auto &&[fidx, s] : c) {
            auto &f = faces()[fidx];
            int region = R[cidx];
            auto &rset = ret[region];
            if (f.is_mesh_face()) {
                rset[s ? 0 : 1].insert(fidx);
            }
        }
    }
    return ret;
}
std::vector<std::array<std::set<int>, 2>> CutCellMesh<3>::orig_face_regions()
    const {
    auto R = face_regions();
    std::vector<std::array<std::set<int>, 2>> ret(R.size());
    for (auto &&[Fsp, OFsp] : mtao::iterator::zip(R, ret)) {
        for (auto &&[Fs, OFs] : mtao::iterator::zip(Fsp, OFsp)) {
            for (auto &&f : Fs) {
                assert(faces()[f].is_mesh_face());
                OFs.insert(faces()[f].as_face_id());
            }
        }
    }
    return ret;
}

mtao::ColVecs3d CutCellMesh<3>::region_centroids() const {
    // delay dividing by 3 until the very end...
    auto R = orig_face_regions();
    mtao::ColVecs3d Cs(3, origF().cols());
    mtao::VecXd vols(origF().cols());
    for (int i = 0; i < Cs.cols(); ++i) {
        auto f = origF().col(i);
        auto c = Cs.col(i);
        mtao::Mat3d A;
        for (int j = 0; j < 3; ++j) {
            A.col(j) = origV().col(f(j));
        }
        vols(i) = A.determinant();
        c = vols(i) * A.rowwise().sum();
    }

    mtao::ColVecs3d C(3, R.size());
    mtao::VecXd V(R.size());
    C.setZero();
    V.setZero();
    for (auto &&[rid, rp] : mtao::iterator::enumerate(R)) {
        auto c = C.col(rid);
        auto &&v = V(rid);
        for (auto &&nf : rp[1]) {
            c -= Cs.col(nf);
            v -= vols(nf);
        }
        for (auto &&pf : rp[0]) {
            c += Cs.col(pf);
            v += vols(pf);
        }
    }

    V = (V.array().abs() > 1e-8).select(V, 1);

    // TODO: handle outside boundaries of the grid using the adaptive grid!
    auto bb = bbox();
    double gv = bb.sizes().prod();
    mtao::Vec3d gc = 1.5 * (bb.min() + bb.min()) * gv;
    C.col(0) += gc;
    V(0) += gv;

    for (int i = 0; i < C.cols(); ++i) {
        C.col(i) /= 3 * V(i);
    }
    return C;
}

CutCellMesh<3> CutCellMesh<3>::from_file(const std::string &filename) {
    return from_proto(filename);
}

void CutCellMesh<3>::serialize(protobuf::CutMeshProto &cmp) const {
    protobuf::serialize(origin(), *cmp.mutable_origin());
    protobuf::serialize(dx(), *cmp.mutable_dx());
    protobuf::serialize(vertex_shape(), *cmp.mutable_shape());

    for (int i = 0; i < cut_vertex_size(); ++i) {
        protobuf::serialize(cut_vertex(i), *cmp.add_vertices());
    }
    // TODO: Add back in after upgrading proto
    /*
        for(int i = 0; i < cut_edge_size(); ++i) {
            protobuf::serialize(m_cut_edges.col(i),*cmp.add_edges());
        }
        */
    for (int i = 0; i < m_origV.cols(); ++i) {
        protobuf::serialize(m_origV.col(i), *cmp.add_origv());
    }
    for (int i = 0; i < m_origF.cols(); ++i) {
        protobuf::serialize(m_origF.col(i), *cmp.add_origf());
    }
    for (auto &&f : m_faces) {
        f.serialize(*cmp.add_faces());
    }
    for (auto &&c : m_cells) {
        c.serialize(*cmp.add_cells());
    }
    for (auto &&i : m_folded_faces) {
        cmp.add_foldedfaces(i);
    }
    auto &&mf = *cmp.mutable_mesh_faces();
    for (auto &&[idx, bmf] : m_mesh_cut_faces) {
        auto &b = mf[idx];
        b.set_parent_id(bmf.parent_fid);
        for (int j = 0; j < bmf.barys.cols(); ++j) {
            protobuf::serialize(bmf.barys.col(j),
                                *b.add_barycentric_coordinates());
        }
    }
    {
        auto &cmap = *cmp.mutable_cubes();
        for (auto &&[c, cell] : m_exterior_grid.cells()) {
            cell.serialize(cmap[c]);
        }
        auto &rmap = *cmp.mutable_cube_regions();
        for (auto &&[a, b] : m_adaptive_grid_regions) {
            rmap[a] = b;
        }
    }
}
CutCellMesh<3> CutCellMesh<3>::from_proto(const std::string &filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (ifs.good()) {
        protobuf::CutMeshProto cmp;
        if (cmp.ParseFromIstream(&ifs)) {
            return from_proto(cmp);
        }
    }
    return {};
}
CutCellMesh<3> CutCellMesh<3>::from_proto(const protobuf::CutMeshProto &cmp) {
    mtao::Vec3d o, dx;
    std::array<int, 3> s;
    o = protobuf::deserialize(cmp.origin());
    dx = protobuf::deserialize(cmp.dx());
    protobuf::deserialize(cmp.shape(), s);

    // s should be the vertex shape
    CutCellMesh<3> ret = CutCellMesh<3>::StaggeredGrid(GridType(s, dx, o));
    // ret.m_face_volumes.resize(20);

    ret.m_cut_vertices.resize(cmp.vertices().size());
    for (int i = 0; i < ret.cut_vertex_size(); ++i) {
        protobuf::deserialize(cmp.vertices(i), ret.m_cut_vertices[i]);
    }
    // TODO: Add back in after updating proto
    /*
        ret.m_cut_edges.resize(2,cmp.edges().size());
        for(int i = 0; i < ret.cut_edge_size(); ++i) {
            ret.m_cut_edges.col(i) = protobuf::deserialize(cmp.edges(i));
        }
        */
    ret.m_origV.resize(3, cmp.origv().size());
    for (int i = 0; i < ret.m_origV.cols(); ++i) {
        ret.m_origV.col(i) = protobuf::deserialize(cmp.origv(i));
    }
    ret.m_origF.resize(3, cmp.origf().size());
    for (int i = 0; i < ret.m_origF.cols(); ++i) {
        ret.m_origF.col(i) = protobuf::deserialize(cmp.origf(i));
    }
    ret.m_faces.resize(cmp.faces().size());
    for (int i = 0; i < cmp.faces().size(); ++i) {
        ret.m_faces[i] = CutFace<3>::from_proto(cmp.faces(i));
        ret.m_faces[i].update_mask(ret.cut_vertices(), ret.vertex_grid());
    }
    ret.m_cells.resize(cmp.cells().size());
    for (int i = 0; i < cmp.cells().size(); ++i) {
        ret.m_cells[i] = CutCell::from_proto(cmp.cells(i));
    }
    std::copy(cmp.foldedfaces().begin(), cmp.foldedfaces().end(),
              std::inserter(ret.m_folded_faces, ret.m_folded_faces.end()));

    for (auto &&[idx, btf] : cmp.mesh_faces()) {
        int pid = btf.parent_id();
        size_t bssize = btf.barycentric_coordinates_size();
        mtao::ColVecs3d B(3, bssize);
        for (int i = 0; i < bssize; ++i) {
            B.col(i) = protobuf::deserialize(btf.barycentric_coordinates(i));
        }
        ret.m_mesh_cut_faces[idx] = {B, pid};
    }

    for (auto &&[a, b] : cmp.cubes()) {
        ret.m_exterior_grid.m_cells[a] = AdaptiveGrid::Cell::from_proto(b);
    }
    ret.m_exterior_grid.make_faces();
    for (auto &&[a, b] : cmp.cube_regions()) {
        ret.m_adaptive_grid_regions[a] = b;
    }

    ret.recompute_active_cell_mask();

    return ret;
}
void CutCellMesh<3>::recompute_active_cell_mask() {
    auto &mask = m_active_grid_cell_mask;
    mask = GridDatab::Constant(false, StaggeredGrid::cell_shape());
    auto C = cell_centroids();
    for (auto &&[i, cell] : mtao::iterator::enumerate(m_cells)) {
        auto cent = C.col(i);
        auto [c, q] = coord(cent);
        mask(c) = true;
    }
}
std::array<mtao::ColVecs2d, 3> CutCellMesh<3>::compute_subVs() const {
    auto V = vertices();
    std::array<mtao::ColVecs2d, 3> subVs;
    for (int d = 0; d < 3; ++d) {
        int n0 = (d + 1) % 3;
        int n1 = (d + 2) % 3;

        auto &&subV = subVs[d];
        subV.resize(2, V.cols());
        subV.row(0) = V.row(n0);
        subV.row(1) = V.row(n1);
    }
    return subVs;
}
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> CutCellMesh<3>::triangulate_face(
    int face_index) const {
    spdlog::warn(
        "Inefficient use of triangulation!  try caching your triangulations");
    if (is_cut_face(face_index)) {
        std::array<mtao::ColVecs2d, 3> subVs = compute_subVs();
        auto [V, F] = m_faces[face_index].triangulate(subVs, true);
        if (V.cols() > 0) {
            return {vertex_grid().world_coord(V), F};
        } else {
            return {mtao::ColVecs3d{}, F};
        }
    } else {
        return {mtao::ColVecs3d{},
                exterior_grid().triangulated_face(face_index,
                                                  exterior_grid_face_offset())};
    }
}

void CutCellMesh<3>::triangulate_faces(bool add_verts) {
    std::array<mtao::ColVecs2d, 3> subVs = compute_subVs();

#if defined(MTAO_TBB_ENABLED)
    tbb::parallel_for_each(
        m_faces.begin(), m_faces.end(),
        [&](CutFace<3> &face) { face.cache_triangulation(subVs, add_verts); });
#else
    int i = 0;
#pragma omp parallel for
    for (i = 0; i < m_faces.size(); ++i) {
        auto &&face = m_faces[i];
        face.cache_triangulation(subVs, add_verts);
    }
#endif
}

mtao::VecXd CutCellMesh<3>::dual_edge_lengths() const {
    return operators::dual_edge_lengths(*this);
}

auto CutCellMesh<3>::face_volumes(bool from_triangulation) const -> VecX {
    return operators::face_volumes(*this, from_triangulation);
}

mtao::VecXd CutCellMesh<3>::mesh_face_mask() const {
    return operators::mesh_face_mask(*this);
}
std::set<int> CutCellMesh<3>::grid_boundary_faces() const {
    return operators::grid_boundary_faces(*this);
}

mtao::VecXd CutCellMesh<3>::primal_hodge2() const {
    return operators::primal_hodge2(*this);
}

mtao::VecXd CutCellMesh<3>::dual_hodge2() const {
    return operators::dual_hodge2(*this);
}
mtao::VecXd CutCellMesh<3>::dual_hodge3() const {
    return operators::dual_hodge3(*this);
}
mtao::VecXd CutCellMesh<3>::primal_hodge3() const {
    return operators::primal_hodge3(*this);
}

Eigen::SparseMatrix<double> CutCellMesh<3>::trilinear_matrix() const {
    return operators::trilinear_matrix(*this);
}
Eigen::SparseMatrix<double> CutCellMesh<3>::face_grid_volume_matrix() const {
    return operators::face_grid_volume_matrix(*this);
}

Eigen::SparseMatrix<double> CutCellMesh<3>::barycentric_matrix() const {
    return operators::barycentric_matrix(*this);
}
Eigen::SparseMatrix<double> CutCellMesh<3>::face_barycentric_volume_matrix()
    const {
    return operators::face_barycentric_volume_matrix(*this);
}
Eigen::SparseMatrix<double> CutCellMesh<3>::boundary(
    bool include_grid_boundary_faces) const {
    return operators::boundary(*this, include_grid_boundary_faces);
}
auto CutCellMesh<3>::active_cell_mask() const -> const GridDatab & {
    return active_grid_cell_mask();
}
auto CutCellMesh<3>::active_cells() const -> std::set<coord_type> {
    auto C = cell_centroids();
    std::set<coord_type> ret;
    for (auto &&[i, cell] : mtao::iterator::enumerate(m_cells)) {
        auto cent = C.col(i);
        auto [c, q] = coord(cent);
        ret.insert(c);
    }
    return ret;
}
size_t CutCellMesh<3>::active_cell_count() const {
    return active_cells().size();
}

std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> CutCellMesh<3>::triangulated_cell(
    int idx, bool use_base, bool use_flap) const {
    if (is_cut_cell(idx)) {
        std::vector<mtao::ColVecs3d> mVs;
        std::vector<mtao::ColVecs3i> mFs;
        for (auto &&[fidx, b] : m_cells[idx]) {
            bool is_flap = m_folded_faces.find(fidx) != m_folded_faces.end();
            bool is_base = !is_flap;
            if ((use_flap && is_flap) || (use_base && !is_flap)) {
                auto &tri = m_faces[fidx].triangulation;
                mtao::ColVecs3d V;
                mtao::ColVecs3i F;
                if (tri) {
                    F = *tri;
                    auto &vert = m_faces[fidx].triangulated_vertices;
                    if (vert) {
                        V = *vert;
                    }
                } else {
                    std::tie(V, F) = triangulate_face(fidx);
                }
                if (!b) {
                    auto tmp = F.row(0).eval();
                    F.row(0) = F.row(1);
                    F.row(1) = tmp;
                }
                mVs.emplace_back(std::move(V));
                mFs.emplace_back(std::move(F));
            }
        }
        return {mtao::eigen::hstack_iter(mVs.begin(), mVs.end()),
                mtao::eigen::hstack_iter(mFs.begin(), mFs.end())};
        /*
            if(mFs.size() > 0) {
                    return
           {mtao::eigen::hstack_iter(mVs.begin(),mVs.end()),mtao::eigen::hstack_iter(mFs.begin(),mFs.end())};
                if(mVs.size() > 0) {
                    return
           {mtao::eigen::hstack_iter(mVs.begin(),mVs.end()),mtao::eigen::hstack_iter(mFs.begin(),mFs.end())};
                } else {
                    return {{},mtao::eigen::hstack_iter(mFs.begin(),mFs.end())};
                }
            }
            */
    } else {
        if (use_base) {
            return {mtao::ColVecs3d{}, m_exterior_grid.triangulated(idx)};
        }
    }
    return {};
}

std::tuple<mtao::ColVecs3d, mtao::ColVecs3i>
CutCellMesh<3>::compact_triangulated_cell(int cell_index) const {
    auto [V, F] = triangulated_cell(cell_index, true, true);
    return mtao::geometry::mesh::compactify(mtao::eigen::hstack(vertices(), V),
                                            F);
}
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i>
CutCellMesh<3>::compact_triangulated_face(int face_index, bool flip) const {
    if (is_cut_face(face_index)) {
        auto &f = m_faces[face_index];
        if (!f.triangulation) {
            mtao::logging::warn() << "Triangulate mesh first!";
            return {};
        } else {
            mtao::ColVecs3d V;
            mtao::ColVecs3i F;
            if (f.triangulated_vertices) {
                std::tie(V, F) = mtao::geometry::mesh::compactify(
                    mtao::eigen::hstack(
                        vertices(),
                        vertex_grid().world_coord(*f.triangulated_vertices)),
                    *f.triangulation);
            } else {
                std::tie(V, F) = mtao::geometry::mesh::compactify(
                    vertices(), *f.triangulation);
            }
            if (flip) {
                auto R = F.row(0).eval();
                F.row(0) = F.row(1);
                F.row(1) = R;
            }
            return {V, F};
        }
    } else {
        auto [V, F] = mtao::geometry::mesh::compactify(
            vertices(), exterior_grid().triangulated_face(
                            face_index, exterior_grid_face_offset()));
        if (flip) {
            auto R = F.row(0).eval();
            F.row(0) = F.row(1);
            F.row(1) = R;
        }
        return {V, F};
    }
}

auto CutCellMesh<3>::cells_by_grid_cell() const
    -> std::map<coord_type, std::set<int>> {
    std::map<coord_type, std::set<int>> R;
    for (auto &&[i, c] : mtao::iterator::enumerate(cells())) {
        R[c.grid_cell].insert(i);
    }
    return R;
}
std::set<int> CutCellMesh<3>::cells_in_grid_cell(const coord_type &c) const {
    std::set<int> R;
    for (auto &&[i, cell] : mtao::iterator::enumerate(cells())) {
        if (cell.grid_cell == c) {
            R.insert(i);
        }
    }
    return R;
}
int CutCellMesh<3>::get_cell_index(const VecCRef &p) const {
    auto v = vertex_grid().local_coord(p);
    // check if its  in an adaptive grid cell, tehn we can just use that cell
    // if the grid returns -2 we pass that through
    if (int ret = m_exterior_grid.get_cell_index(v); ret != -1) {
        if (ret == -2) {
            mtao::logging::warn() << "Point lies outside the grid";
        }
        return ret;
    } else {
        auto [c, q] = vertex_grid().coord(p);
        auto cell_indices = cells_in_grid_cell(c);
        ColVecs V;
        const auto &CV = cached_vertices();
        const ColVecs *VP = CV ? &*CV : &V;
        if (!CV) {
            V = vertices();
        }
        for (auto &&ci : cell_indices) {
            auto &&cell = cells().at(ci);
            if (cell.contains(*VP, faces(), p)) {
                return ci;
            }
        }
        // TODO
        // if I haven't returned yet then either this ccm is bad
        // or i'm too close to an edge. lets not assume that for now
        mtao::logging::warn()
            << "There are CCM in this grid cell but I failed to find it!"
            << c[0] << " " << c[1] << " " << c[2];
    }
    return -1;
}

bool CutCellMesh<3>::is_in_cell(const VecCRef &p, size_t index) const {
    auto v = vertex_grid().local_coord(p);
    // check if its  in an adaptive grid cell, tehn we can just use that cell
    // if the grid returns -2 we pass that through
    if (int ret = m_exterior_grid.get_cell_index(v); ret != -1) {
        if (ret == -2) {
            mtao::logging::warn() << "Point lies outside the grid";
        }
        return ret == index;
    } else {
        auto &&cell = cells().at(index);
        ColVecs V;
        const auto &CV = cached_vertices();
        const ColVecs *VP = CV ? &*CV : &V;
        if (!CV) {
            V = vertices();
        }
        if (cell.contains(*VP, faces(), p)) {
            return true;
        }
    }
    // TODO
    // if I haven't returned yet then either this ccm is bad
    // or i'm too close to an edge. lets not assume that for now
    return false;
}

auto CutCellMesh<3>::edges() const -> Edges {
    // TODO: this is not a real set of edges. Gotta figure out what i really
    // want to do here

    Edges VE =
        mtao::geometry::grid::GridTriangulator<GridType>(vertex_grid()).edges();
    Edges EE = cut_edges_eigen();
    return mtao::eigen::hstack(VE, EE);
    // return Base::edges();
}
}  // namespace mandoline
