#include <igl/AABB.h>

#include "mandoline/mesh3.hpp"
#include "mandoline/operators/boundary3.hpp"
#include "mandoline/operators/interpolation3.hpp"
#include "mandoline/operators/cell_indices.hpp"
#include "mandoline/operators/nearest_facet.hpp"
#include "mtao/geometry/winding_number.hpp"
namespace mandoline::operators {
namespace {
coord_mask<3> projection_mask(const auto& grid, const auto& cell) {
    coord_mask<3> mask;
    auto s = grid.shape();
    for (auto&& [b, c, m] : mtao::iterator::zip(s, cell, mask)) {
        if (c < 0) {
            m = 0;
        }
        if (c >= b) {
            m = b - 1;
        }
    }
    return mask;
}
}  // namespace

BoundaryFacetProjector3::BoundaryFacetProjector3(const CutCellMesh<3>& ccm)
    : Base(static_cast<const Base&>(ccm)) {
    _V = ccm.origV().transpose();
    _F = ccm.origF().transpose();
    _aabb = std::make_shared<igl::AABB<mtao::RowVecs3d, 3>>();
    _aabb->init(_V, _F);
    _cell_ownership_grid = ccm.exterior_grid().cell_ownership_grid();
}
BoundaryFacetProjector3::~BoundaryFacetProjector3() {}

std::tuple<int, double> BoundaryFacetProjector3::nearest_triangle(
    const Eigen::Ref<const mtao::Vec3d> p) const {
    std::tuple<int, double> ret;
    auto& [idx, dist] = ret;
    mtao::RowVector<double, 3> c;
    dist = std::sqrt(_aabb->squared_distance(_V, _F, p.transpose(), idx, c));
    return ret;
}
std::tuple<int, double> BoundaryFacetProjector3::nearest_grid_face(
    const Eigen::Ref<const mtao::Vec3d> p) const {
    std::tuple<int, double> ret = {-1, std::numeric_limits<double>::max()};
    auto& [idx, dist] = ret;
    int axis = 0;
    using coord_type = std::array<int, 3>;

    auto [grid_cell, quotient] = vertex_grid().coord(p);
    auto mask = projection_mask(cell_grid(), grid_cell);
    const int facet_case = mask.count();

    coord_type projected_grid_cell = grid_cell;
    mask.clamp(projected_grid_cell);
    // std::cout << facet_case << "(" << std::string(mask) << std::endl;

    /*
    idx = facet_case;

        auto bb = bbox();
        if(bb.contains(p)) {
            dist = 88888;
        } else {
        auto dM = (p - bb.max());
        auto dm = (bb.min() - p);
        mtao::Vec3d diff = (dM.array() > 0).select(dM, 0);
        diff = (dm.array() > 0).select(dm, diff);
        dist = diff.norm();
        }
        */
    if (facet_case == 0) {  // in a grid cell
        bool above = true;
        for (int j = 0; j < 3; ++j) {
            if (quotient[j] < dist) {
                dist = quotient[j];
                axis = j;
                above = false;
            }
            double tv = 1 - quotient[j];
            if (tv < dist) {
                dist = tv;
                axis = j;
                above = true;
            }
        }
        if (above) {
            projected_grid_cell[axis] += 1;
        }
        dist /= dx()(axis);
    } else if (facet_case == 1) {  // in a grid face
        axis = mask.bound_axis();
        // spdlog::info("Doing face case on axis {}", axis);
        dist = quotient[axis] + grid_cell[axis] - projected_grid_cell[axis];
        if (dist >= 0) {
            projected_grid_cell[axis] += 1;
        }
        dist = std::abs(dist);
        dist /= dx()(axis);
    } else if (facet_case == 2) {
        int edge_axis = mask.unbound_axis();
        // spdlog::info("Edge axis: {} with pc {}", edge_axis,
        //             fmt::join(projected_grid_cell, ","));
        bool above;
        mtao::Vec2d local_dists;
        for (int j = 0; j < 2; ++j) {
            int my_axis = (j + edge_axis + 1) % 3;
            double& my_dist = local_dists(j) = quotient[my_axis] +
                                               grid_cell[my_axis] -
                                               projected_grid_cell[my_axis];
            double adist = std::abs(my_dist);
            if (adist < dist) {
                dist = adist;
                axis = my_axis;
                // spdlog::info("axis {} got dist {}", dist, axis);
                above = my_dist >= 1;
            }
            my_dist /= dx()(my_axis);
        }
        if (above) {
            projected_grid_cell[axis] += 1;
        }
        dist = local_dists.norm();
    }
    if (facet_case == 3) {
        bool above;
        for (int j = 0; j < 3; ++j) {
            double my_dist =
                quotient[j] + grid_cell[j] - projected_grid_cell[j];
            double adist = std::abs(my_dist);
            if (adist < dist) {
                dist = adist;
                axis = j;
                above = my_dist >= 0;
            }
        }
        if (above) {
            projected_grid_cell[axis] += 1;
        }

        auto bb = bbox();
        auto dM = (p - bb.max());
        auto dm = (bb.min() - p);
        mtao::Vec3d diff = (dM.array() > 0).select(dM, 0);
        diff = (dm.array() > 0).select(dm, diff);
        dist = diff.norm();
    }
    idx = staggered_index<2>(projected_grid_cell, axis);
    // spdlog::info("Staggered index of {} on axis {} => {}",
    //             fmt::join(projected_grid_cell, ","), axis, idx);
    return ret;
}

std::vector<std::tuple<int, double>> BoundaryFacetProjector3::nearest_triangles(
    const Eigen::Ref<const mtao::ColVecs3d> P) const {
    std::vector<std::tuple<int, double>> I(P.cols());
    tbb::parallel_for<int>(0, I.size(), [&](int idx) {
        auto& pr = I[idx];
        auto p = P.col(idx);
        pr = nearest_triangle(p);
    });
    return I;
}
std::vector<std::tuple<int, double>>
BoundaryFacetProjector3::nearest_grid_faces(
    const Eigen::Ref<const mtao::ColVecs3d> P) const {
    std::vector<std::tuple<int, double>> I(P.cols());
    tbb::parallel_for<int>(0, I.size(), [&](int idx) {
        auto& pr = I[idx];
        auto p = P.col(idx);
        pr = nearest_grid_face(p);
    });
    return I;
}

CellParentMaps3::~CellParentMaps3() {}
CellParentMaps3::CellParentMaps3(const CutCellMesh<3>& ccm) : _projector(ccm) {
    grid_contained_cells.resize(ccm.cell_grid().size());
    for (auto&& [cell_index, cell] : mtao::iterator::enumerate(ccm.cells())) {
        grid_contained_cells[ccm.StaggeredGrid::cell_index(cell.grid_cell)]
            .emplace(cell_index);
    }
#if defined(MANDOLINE_USE_ADAPTIVE_GRID)
    for (auto&& [cid, desc] : ccm.exterior_grid().cells()) {
        auto&& abc = desc.corner();
        auto&& jump = desc.width();
        for (int aa = abc[0]; aa < abc[0] + jump; ++aa) {
            for (int bb = abc[1]; bb < abc[1] + jump; ++bb) {
                for (int cc = abc[2]; cc < abc[2] + jump; ++cc) {
                    int grid_index = ccm.StaggeredGrid::cell_index(
                        std::array<int, 3>{{aa, bb, cc}});
                    grid_contained_cells[grid_index].emplace(cid);
                }
            }
        }
    }
#else
    for (auto&& [cid, coord] :
         mtao::iterator::enumerate(ccm.exterior_grid().cell_coords())) {
        int grid_index = ccm.StaggeredGrid::cell_index(coord);
        grid_contained_cells[grid_index].emplace(cid);
    }

#endif

    cut_cell_coboundary.resize(ccm.face_size());
    {
        auto B = boundary(ccm, true);
        for (int o = 0; o < B.outerSize(); ++o) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(B, o); it;
                 ++it) {
                cut_cell_coboundary[it.row()].emplace(it.col());
            }
        }
    }
    /*
    {
        int count = 0;
        for (auto&& f : grid_contained_faces) {
            if (f.size() == 0) {
                count++;
            }
        }
        spdlog::info("Cut-faces filled ind {} of {} grid faces", count,
                     grid_contained_faces.size());
    }
    */

    grid_contained_faces.resize(ccm.form_size<2>());
    {
        auto G = mandoline::operators::face_grid_volume_matrix(ccm, true);
        for (int o = 0; o < G.outerSize(); ++o) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(G, o); it;
                 ++it) {
                grid_contained_faces[it.col()].emplace(it.row());
            }
        }
    }

    triangle_contained_faces.resize(ccm.origF().cols());
    for (auto&& [face_index, face] :
         mtao::iterator::enumerate(ccm.cut_faces())) {
        if (face.is_mesh_face()) {
            triangle_contained_faces[face.as_face_id()].emplace(face_index);
        } else {
            continue;
            const auto& cells = cut_cell_coboundary[face_index];

            std::array<int, 3> coord{{0, 0, 0}};
            for (auto&& c : cells) {
                // should only have two boundary grid cells and we want the one
                // with higher coordinates
                if (ccm.is_cut_cell(c)) {
                    const auto& gc = ccm.cells()[c].grid_cell;
                    // spdlog::info("Got grid cell {} from index {} / {}",
                    //             fmt::join(gc, ","), c, ccm.cell_size());
                    for (auto&& [a, b] : mtao::iterator::zip(coord, gc)) {
                        a = std::max(a, b);
                    }
                } else {
                    // exterior grid faces
                }
            }
            // spdlog::info("Got resulting cell {} / {}", fmt::join(coord, ","),
            //             fmt::join(ccm.cell_shape(), ","));
            int index = ccm.staggered_index<2>(coord, face.as_axial_axis());
            // spdlog::info("Resulting index {} / {} <- {}", index,
            //             grid_contained_faces.size(), face_index);
            grid_contained_faces[index].emplace(face_index);
        }
    }
    /*
    {
        int count = 0;
        for (auto&& f : grid_contained_faces) {
            if (f.size() == 0) {
                count++;
            }
        }
        spdlog::info("From just cut-faces i filled {} of {} grid faces", count,
                     grid_contained_faces.size());
    }
    int cf_offset = ccm.cut_face_size();
    for (auto&& [fidx_, face] :
         mtao::iterator::enumerate(ccm.exterior_grid().faces())) {
        const int fidx = fidx_ + cf_offset;
        const int width = face.width();
        const int axis = face.axis();
        const std::array<int, 3>& corner = face.corner();

        const int& cornu = corner[(axis + 1) % 3];
        const int& cornv = corner[(axis + 2) % 3];
        std::array<int, 3> c = corner;
        int& cu = c[(axis + 1) % 3];
        int& cv = c[(axis + 2) % 3];
        for (cu = cornu; cu < cornu + 1; ++cu) {
            for (int cv = cornv; cv < cornv + 1; ++cv) {
                grid_contained_faces[ccm.StaggeredGrid::staggered_index<2>(
                                         c, axis)]
                    .emplace(fidx);
            }
        }
    }
    for (auto&& it :
         ccm.exterior_grid().grid_face_projection(ccm.cut_face_size())) {
        grid_contained_faces[it.col()].emplace(it.row());
        spdlog::info("{} {}", it.row(), it.col());
    }
    */

    subVs = ccm.compute_subVs();
}

/*
namespace {

int nearest_cut_vertex(const CutCellMesh<3> &ccm,
                   Eigen::Ref<const mtao::Vec3d> p,
                   const coord_mask<3> &grid_case) {}
int nearest_cut_edge(const CutCellMesh<3> &ccm, Eigen::Ref<const mtao::Vec3d> p,
                 const coord_mask<3> &grid_case) {
int edge_index = -1;
std::array<int, 2> edge_endpoints;
for (auto &&[eidx, edge] : mtao::iterator::enumerate(cut_edges())) {
    if (edge.is_axial_edge()) {
        int axis = edge.as_axial_axis();
        if (edge.mask() == mask &&
            edge.as_axial_coord() == grid_cell[axis]) {
            double t = quotient[axis];
            auto [a, b] = edge.indices;
            if (t >= masked_vertex(a).quot(axis) &&
                t < masked_vertex(b).quot(axis)) {
                edge_index = eidx;
                edge_endpoints = edge.indices;
            }
        }
    }
}
if (edge_index < 0) {
    throw std::runtime_error(
        "Cutcell mesh wasnt made properly because we couldnt find "
        "an edge");
}
}

int nearest_cut_face(const CutCellMesh<3> &ccm,
                 Eigen::Ref<const mtao::Vec3d> &p,
                 const coord_mask<3> &grid_case) {
for (auto &&cell_index : cell_indices) {
    for (auto &&[fidx, i] : cells().at(cell_index)) {
        const auto &face = cut_face(fidx);
        if (bool(face.external_boundary) &&
            std::get<0>(*face.external_boundary) == -2) {
            if (face.indices.size() != 1) {
                throw std::runtime_error(
                    "projecting to nontrivial boundary faces not "
                    "implemented yet");
            }
            for (auto &&v : *face.indices.begin()) {
                if (!is_grid_vertex(v)) {
                    throw std::runtime_error(
                        "projecting non-gridlike faces not "
                        "implemented yet");
                }
            }
            if (face.mask() == mask) {
                return cell_index;
            }
        }
    }
}
}

int nearest_cut_cell(const CutCellMesh<3> &ccm, Eigen::Ref<const mtao::Vec3d> p,
                 const coord_mask<3> &grid_case) {
int count = grid_case.count();
coord_type projected_grid_cell = grid_cell;
mask.clamp(projected_grid_cell);

auto cell_indices = cut_cells_in_grid_cell(grid_cell);
if (count == 3) {
    if (cell_indices.size() == 1) {
        return *cell_indices.begin();
    } else {
        spdlog::warn(
            "Get nearest cell index failed because mtao was too "
            "lazy to implement the vertex case");
        return -1;
    }
} else if (count == 2) {
    int edge_index = nearest_edge(ccm, p);
    for (auto &&cell_index : cell_indices) {
        for (auto &&[fidx, i] : cells().at(cell_index)) {
            const auto &face = cut_face(fidx);
            if (bool(face.external_boundary) &&
                std::get<0>(*face.external_boundary) == -2) {
                for (auto &&loop : face.indices) {
                    for (size_t j = 0; j < loop.size(); ++j) {
                        const auto &[a, b] = edge_endpoints;
                        const auto &c = loop[j];
                        const auto &d = loop[(j + 1) % loop.size()];
                        if ((a == c && b == d) || (a == d && b == c)) {
                            return cell_index;
                        }
                    }
                }
            }
        }
    }
    if (edge_index < 0) {
        throw std::runtime_error(fmt::format(
            "Cutcell mesh wasnt made properly because we couldnt find "
            "a face that uses cutedge {}",
            edge_index));
    }

} else if (count == 1) {
    int face_index = nearest_face(ccm, p);
} else if (count == 0) {
    spdlog::error(
        "get_nearest_cell_index can only get a cell interior during "
        "projection if the mesh was not produced properly!");
}
}
else {
return index;
}
return -1;
}  // namespace

mtao::VecXi nearest_vertices(const CutCellMesh<3> &ccm,
                         Eigen::Ref<const mtao::ColVecs3d> p);
mtao::VecXi nearest_edges(const CutCellMesh<3> &ccm,
                      Eigen::Ref<const mtao::ColVecs3d> p);
mtao::VecXi nearest_faces(const CutCellMesh<3> &ccm,
                      Eigen::Ref<const mtao::ColVecs3d> p);
*/
mtao::VecXi nearest_faces(const CutCellMesh<3>& ccm,
                          Eigen::Ref<const mtao::ColVecs3d> P) {
    CellParentMaps3 parent_maps(ccm);
    return nearest_faces(ccm, parent_maps, P);
}

mtao::VecXi nearest_faces(const CutCellMesh<3>& ccm,
                          const CellParentMaps3& parent_maps,
                          Eigen::Ref<const mtao::ColVecs3d> P) {
    auto nearest_triangles = parent_maps._projector.nearest_triangles(P);
    auto nearest_grid_faces = parent_maps._projector.nearest_grid_faces(P);
    mtao::VecXi I(P.cols());
    I.setConstant(-1);

    const auto& subVs = parent_maps.subVs;
    // tbb::parallel_for<int>(0, I.size(), [&](int j) {
    for (int j = 0; j < I.size(); ++j) {
        const auto& [parent_triangle, tri_distance] = nearest_triangles[j];
        const auto& [parent_grid_face, grid_distance] = nearest_grid_faces[j];

        auto p = P.col(j);
        auto& index = I(j);

        std::cout << tri_distance << " " << grid_distance << std::endl;
        if (false && tri_distance < grid_distance) {
            auto tf = ccm.origF().col(parent_triangle);
            auto a = ccm.origV().col(tf(0));
            auto b = ccm.origV().col(tf(1));
            auto c = ccm.origV().col(tf(2));
            mtao::Matrix<double, 3, 2> B;
            B.col(0) = b - a;
            B.col(1) = c - a;
            // A [1,0] = B
            //   [0,1]
            // B^T A = B^T B
            // (B^T B )^{-1} B^T A = I
            // (B^T B)^{-1} B^T = A^{-1}
            // A x = P
            // x = (B^T B)^{-1} B^T P

            mtao::Vec2d p2 =
                (B.transpose() * B).inverse() * (B.transpose() * p);
            bool move_x = p2.x() <= 0;
            bool move_y = p2.y() <= 0;
            bool move_z = p2.sum() >= 1;
            if (move_x && move_y) {
                p2.setZero();
            } else if (move_x && move_z) {
                p2.setUnit(1);
            } else if (move_y && move_z) {
                p2.setUnit(0);
            } else if (move_x) {
                p2.y() /= 1 - p2.x();  // p2.y() / (1 - p2.sum()) + p2.y()) =
                                       // p2.y() / (1-p2.x());
                p2.x() = 0;
            } else if (move_y) {
                p2.x() /= 1 - p2.y();
            } else if (move_z) {
                p2 /= p2.sum();
            }

            // spdlog::info("Triangle contained faces got {} / {}",
            // parent_triangle, parent_maps.triangle_contained_faces.size());
            const auto& cutfaces =
                parent_maps.triangle_contained_faces[parent_triangle];
            if (cutfaces.size() == 1) {
                index = *cutfaces.begin();
            } else if (cutfaces.size() > 1) {
                for (auto&& cutface_index : cutfaces) {
                    const auto& mesh_face =
                        ccm.mesh_cut_faces().at(cutface_index);
                    std::vector<int> loop(mesh_face.barys.cols());
                    std::iota(loop.begin(), loop.end(), 0);

                    if (mtao::geometry::interior_winding_number(
                            mesh_face.barys.bottomRows<2>(), loop, p2)) {
                        index = cutface_index;
                        break;
                    }
                }
                if (index == -1) {  // almost definitely a boundary value
                    for (auto&& cutface_index : cutfaces) {
                        const auto& mesh_face =
                            ccm.mesh_cut_faces().at(cutface_index);
                        const auto& B = mesh_face.barys;
                        for (int j = 0; j < B.cols(); ++j) {
                            auto ba = B.col(j);
                            auto bb = B.col((j + 1) % B.cols());
                            double tval;
                            double a, b;
                            if (move_x) {
                                if (ba.x() == 0 && bb.x() == 0) {
                                    tval = p2.y();
                                    std::tie(a, b) =
                                        std::make_tuple(ba.y(), bb.y());
                                } else {
                                    continue;
                                }
                            } else if (move_y) {
                                if (ba.y() == 0 && bb.y() == 0) {
                                    tval = p2.x();
                                    std::tie(a, b) =
                                        std::make_tuple(ba.x(), bb.x());
                                } else {
                                    continue;
                                }
                            } else if (move_z) {
                                if (ba.z() == 0 && bb.z() == 0) {
                                    tval = p2.x();
                                    std::tie(a, b) =
                                        std::make_tuple(ba.x(), bb.x());
                                } else {
                                    continue;
                                }
                            }
                            if (b < a) {
                                std::swap(a, b);
                            }
                            if (a <= tval && tval <= b) {
                                index = cutface_index;
                                break;
                            }
                        }
                        if (index >= 0) {
                            break;
                        }
                    }
                }
                if (index == -1) {
                    spdlog::error(
                        "Failed to find a mesh cutface holding {} in triangle "
                        "{} in nearest_faces",
                        fmt::join(p, ","), parent_triangle);
                }
            } else {
                spdlog::error(
                    "Failed to find a mesh cutface holding {} in triangle {} "
                    "due to 0 faces",
                    fmt::join(p, ","), parent_triangle);
            }
        } else {
            const auto& nearest_cut_faces =
                parent_maps.grid_contained_faces[parent_grid_face];
            if (nearest_cut_faces.size() == 1) {
                int cf = *nearest_cut_faces.begin();
                index = cf;
                continue;  // return;
            } else {
                int axis = ccm.form_type<2>(parent_grid_face);
                mtao::Vec2d pp;
                pp(0) = p((axis + 1) % 3);
                pp(1) = p((axis + 2) % 3);
                const auto& V2 = subVs[axis];
                auto is_inside = [&](auto&& face) -> bool {
                    double wn = 0;
                    for (auto&& c : face.indices) {
                        double mywn = mtao::geometry::winding_number(V2, c, pp);
                        wn += mywn;
                    }
                    return std::abs(wn) > .5;
                };
                const auto& nearest_grid_face_cut_faces =
                    parent_maps.grid_contained_faces[parent_grid_face];
                for (auto&& fidx : nearest_grid_face_cut_faces) {
                    if (is_inside(ccm.cut_face(fidx))) {
                        index = fidx;
                        // return;
                        continue;  // return;
                    }
                }

                if (index == -1) {
                    spdlog::error(
                        "Failed to find a grid cutface holding {} in face {}",
                        fmt::join(p, ","), parent_grid_face);
                }
            }
        }
        //});
    }
    return I;
}

mtao::VecXi nearest_mesh_cut_faces(const CutCellMesh<3>& ccm,
                                   Eigen::Ref<const mtao::ColVecs3d> P) {
    CellParentMaps3 parent_maps(ccm);
    return nearest_mesh_cut_faces(ccm, parent_maps, P);
}
mtao::VecXi nearest_mesh_cut_faces(const CutCellMesh<3>& ccm,
                                   const CellParentMaps3& parent_maps,
                                   Eigen::Ref<const mtao::ColVecs3d> P) {
    auto nearest_triangles = parent_maps._projector.nearest_triangles(P);
    mtao::VecXi I(P.cols());
    I.setConstant(-1);
    tbb::parallel_for<int>(0, I.size(), [&](int j) {
        const auto& [parent_triangle, distance] = nearest_triangles[j];
        auto p = P.col(j);
        auto& index = I(j);

        auto tf = ccm.origF().col(parent_triangle);
        auto a = ccm.origV().col(tf(0));
        auto b = ccm.origV().col(tf(1));
        auto c = ccm.origV().col(tf(2));
        mtao::Matrix<double, 3, 2> B;
        B.col(0) = b - a;
        B.col(1) = c - a;
        // A [1,0] = B
        //   [0,1]
        // B^T A = B^T B
        // (B^T B )^{-1} B^T A = I
        // (B^T B)^{-1} B^T = A^{-1}
        // A x = P
        // x = (B^T B)^{-1} B^T P

        mtao::Vec2d p2 = (B.transpose() * B).inverse() * (B.transpose() * p);
        bool move_x = p2.x() <= 0;
        bool move_y = p2.y() <= 0;
        bool move_z = p2.sum() >= 1;
        if (move_x && move_y) {
            p2.setZero();
        } else if (move_x && move_z) {
            p2.setUnit(1);
        } else if (move_y && move_z) {
            p2.setUnit(0);
        } else if (move_x) {
            p2.y() /= 1 - p2.x();  // p2.y() / (1 - p2.sum()) + p2.y()) = p2.y()
                                   // / (1-p2.x());
            p2.x() = 0;
        } else if (move_y) {
            p2.x() /= 1 - p2.y();
        } else if (move_z) {
            p2 /= p2.sum();
        }

        // spdlog::info("Triangle contained faces got {} / {}", parent_triangle,
        // parent_maps.triangle_contained_faces.size());
        const auto& cutfaces =
            parent_maps.triangle_contained_faces[parent_triangle];
        if (cutfaces.size() == 1) {
            index = *cutfaces.begin();
        } else if (cutfaces.size() > 1) {
            for (auto&& cutface_index : cutfaces) {
                const auto& mesh_face = ccm.mesh_cut_faces().at(cutface_index);
                std::vector<int> loop(mesh_face.barys.cols());
                std::iota(loop.begin(), loop.end(), 0);

                if (mtao::geometry::interior_winding_number(
                        mesh_face.barys.bottomRows<2>(), loop, p2)) {
                    index = cutface_index;
                    break;
                }
            }
            if (index == -1) {  // almost definitely a boundary value
                for (auto&& cutface_index : cutfaces) {
                    const auto& mesh_face =
                        ccm.mesh_cut_faces().at(cutface_index);
                    const auto& B = mesh_face.barys;
                    for (int j = 0; j < B.cols(); ++j) {
                        auto ba = B.col(j);
                        auto bb = B.col((j + 1) % B.cols());
                        double tval;
                        double a, b;
                        if (move_x) {
                            if (ba.x() == 0 && bb.x() == 0) {
                                tval = p2.y();
                                std::tie(a, b) =
                                    std::make_tuple(ba.y(), bb.y());
                            } else {
                                continue;
                            }
                        } else if (move_y) {
                            if (ba.y() == 0 && bb.y() == 0) {
                                tval = p2.x();
                                std::tie(a, b) =
                                    std::make_tuple(ba.x(), bb.x());
                            } else {
                                continue;
                            }
                        } else if (move_z) {
                            if (ba.z() == 0 && bb.z() == 0) {
                                tval = p2.x();
                                std::tie(a, b) =
                                    std::make_tuple(ba.x(), bb.x());
                                // spdlog::info("Intervaling {} between
                                // [{},{}]",
                                //             tval, a, b);
                            } else {
                                continue;
                            }
                        }
                        if (b < a) {
                            std::swap(a, b);
                        }
                        if (a <= tval && tval <= b) {
                            index = cutface_index;
                            break;
                        }
                    }
                    if (index >= 0) {
                        break;
                    }
                }
                if (index == -1) {
                    spdlog::error(
                        "Failed to find a mesh cutface holding {} in triangle "
                        "{} through edges (bary was {})",
                        fmt::join(p, ","), parent_triangle, fmt::join(p2, ","));
                }
            }
        } else {
            spdlog::error(
                "Failed to find a mesh cutface holding {} in triangle {}",
                fmt::join(p, ","), parent_triangle);
        }
    });
    return I;
}

mtao::VecXi nearest_cells(const CutCellMesh<3>& ccm,
                          Eigen::Ref<const mtao::ColVecs3d> P) {
    CellParentMaps3 parent_maps(ccm);
    return nearest_cells(ccm, parent_maps, P);
}
mtao::VecXi nearest_cells(const CutCellMesh<3>& ccm,
                          const CellParentMaps3& parent_maps,
                          Eigen::Ref<const mtao::ColVecs3d> P) {
    mtao::VecXi I(P.cols());
    // I.setConstant(0);
    // return I;

    auto nearest_grid_faces = parent_maps._projector.nearest_grid_faces(P);
    // auto nearest_triangle_faces =
    // parent_maps._projector.nearest_triangle_faces(P);
    auto bbox = ccm.bbox();

    // auto is_boundary_grid_face = [&](int index) {
    //    int axis = ccm.form_type<2>(index);
    //    coord_type coord =
    //    ccm.staggered_unindex<2>(nearest_grid_index,axis); int val =
    //    coord[axis]; return val == 0 || coord = ccm.shape()[axis];

    //}:
    const auto& subVs = parent_maps.subVs;

    // vertices rewritten as projected to each cell
    mtao::ColVecs3d _V;
    const auto& CV = ccm.cached_vertices();
    const mtao::ColVecs3d* V = CV ? &*CV : &_V;
    if (!CV) {
        _V = ccm.vertices();
    }
#if defined(MTAO_TBB_ENABLED)
    tbb::parallel_for<int>(0, I.size(), [&](const int j) {
#else
    for (int j = 0; j < I.size(); ++j) {
#endif
        auto p = P.col(j);
        int& index = I(j);
        auto [grid_cell, quotient] = ccm.vertex_grid().coord(p);
        // if we're on the inside then use the gotten cell
        if (bbox.contains(p)) {
            int got_cell = get_cell_index(ccm,parent_maps._projector._cell_ownership_grid,p);
            if (got_cell >= 0) {
                index = got_cell;
#if defined(MTAO_TBB_ENABLED)
            return;
#else
            continue;
#endif
            }
        }

        const auto [nearest_grid_index, grid_distance] = nearest_grid_faces[j];
        const auto& nearest_cut_faces =
            parent_maps.grid_contained_faces[nearest_grid_index];
        // spdlog::info("Nearest grid face {} has {} cut-faces",
        //             nearest_grid_index, nearest_cut_faces.size());
        if (nearest_cut_faces.size() == 1) {
            int cf = *nearest_cut_faces.begin();
            const auto& cells = parent_maps.cut_cell_coboundary[cf];
            // spdlog::info("Cut-face {} has {} coboundaries", cf,
            // cells.size());
            if (cells.size() == 1) {
                index = *cells.begin();
#if defined(MTAO_TBB_ENABLED)
            return;
#else
            continue;
#endif
            }
        }
        std::cout << std::endl;
        spdlog::warn(
            "Exterior grid cases! vertex {} has local coordinate {} + {}",
            fmt::join(p, ","), fmt::join(grid_cell, ","),
            fmt::join(quotient, ","));
        parent_maps._projector.nearest_grid_face(p);
        // const auto [nearest_triangle_index,triangle_distance] =
        // nearest_triangle_faces[j];
        const int axis = ccm.form_type<2>(nearest_grid_index);
        std::array<int, 3> coord =
            ccm.staggered_unindex<2>(nearest_grid_index, axis);
        const bool near_boundary = coord[axis] == 0;
        const bool far_boundary = coord[axis] == ccm.cell_shape()[axis];
        spdlog::info("Case near {} far {} from coord {} / {} on axis {}",
                     near_boundary, far_boundary, fmt::join(coord, ","),
                     fmt::join(ccm.cell_shape(), ","), axis);
        const std::set<int>* potential_cells_ptr;
        if (near_boundary) {
            potential_cells_ptr =
                &parent_maps.grid_contained_cells[ccm.StaggeredGrid::cell_index(
                    coord)];
            spdlog::info("USing cells from grid cell {}, which has {}",
                         fmt::join(coord, ","),
                         fmt::join(*potential_cells_ptr, ","));
        } else if (far_boundary) {
            coord[axis]--;
            spdlog::info("USing cells from grid cell {}, which has {}",
                         fmt::join(coord, ","),
                         fmt::join(*potential_cells_ptr, ","));
            potential_cells_ptr =
                &parent_maps.grid_contained_cells[ccm.StaggeredGrid::cell_index(
                    coord)];
        } else {
            // auto [grid_cell, quotient] = ccm.vertex_grid().coord(p);
            // for (auto&& [i, s] :
            //        mtao::iterator::zip(grid_cell,
            //        ccm.cell_grid().shape())) {
            //    i = std::clamp(i, 0, s);
            //    potential_cells_ptr =
            //        &parent_maps.grid_contained_cells
            //        [ccm.StaggeredGrid::cell_index(coord)];
            //}
            spdlog::error(
                "Interior faces should find neighbors with get_cell; "
                "setting "
                "particle {} to cell 0",
                j);
            index = 0;
#if defined(MTAO_TBB_ENABLED)
            return;
#else
            continue;
#endif
        }
        const std::set<int>& potential_cells = *potential_cells_ptr;

        if (potential_cells.size() == 0) {
            spdlog::error("No possible cells! giving index 0");
            index = 0;
#if defined(MTAO_TBB_ENABLED)
            return;
#else
            continue;
#endif
        } else if (potential_cells.size() == 1) {
            index = *potential_cells.begin();
#if defined(MTAO_TBB_ENABLED)
            return;
#else
            continue;
#endif
        } else {
            auto [grid_cell, quotient] = ccm.vertex_grid().coord(p);
            auto mask = projection_mask(ccm.cell_grid(), grid_cell);
            int facet_case = mask.count();
            index = -1;

            // try processing faces
            if (facet_case == 1) {
                int axis = mask.bound_axis();
                mtao::Vec2d pp;
                pp(0) = p((axis + 1) % 3);
                pp(1) = p((axis + 2) % 3);
                const auto& V2 = subVs[axis];
                auto is_inside = [&](auto&& face) -> bool {
                    double wn = 0;
                    for (auto&& c : face.indices) {
                        double mywn = mtao::geometry::winding_number(V2, c, pp);
                        wn += mywn;
                    }
                    return std::abs(wn) > .5;
                };
                const auto& nearest_grid_face_cut_faces =
                    parent_maps.grid_contained_faces[nearest_grid_index];
                for (auto&& fidx : nearest_grid_face_cut_faces) {
                    if (is_inside(ccm.cut_face(fidx))) {
                        const auto& cells =
                            parent_maps.cut_cell_coboundary.at(fidx);
                        if (cells.size() == 1) {
                            index = *cells.begin();
                        } else {
                            spdlog::warn(
                                "Boundary faces cant have multiple "
                                "interior "
                                "({} I see )"
                                "cells. picking the first one",
                                cells.size());
                            index = *cells.begin();
                        }
                        break;
                    }
                }
                if (index != -1) {
                    spdlog::error("Failed to find nearest grid cutface");
                }
            }
            // spdlog::warn(
            //    "Did not implement non-trivial edge and vertex cell
            //    projection " "cases! vertex {} has local coordinate {} +
            //    {} with {} options "
            //    "({})",
            //    fmt::join(p, ","), fmt::join(coord, ","),
            //    fmt::join(quotient, ","), potential_cells.size(),
            //    fmt::join(potential_cells, ","));
            index = *potential_cells.begin();

            /*
            // let stragglers from the face case come along
            if (facet_case <= 2) {
                auto point_edge_distance = [&](auto&& cutedge) -> bool {
                    auto&& [ai, bi] = cutedge.indices;
                    auto a = V.col(ai);
                    auto b = V.col(bi);
                    auto T = (b - a).eval();
                    double D = T.norm();
                    double term = std::clamp<double>((T.dot(p - a) / D), 0,
            1); auto pp = a + term * T; return (p - pp).norm();
                };

            }
            */
        }
#if defined(MTAO_TBB_ENABLED)
    });
#else
    }
#endif
    return I;
    /*

    const auto g = exterior_grid().cell_ownership_grid();
    const auto grid_to_cells = cut_cells_per_grid_cell();
    ColVecs V;
    const auto& CV = cached_vertices();
    const ColVecs* VP = CV ? &*CV : &V;
    if (!CV) {
        V = vertices();
    }
//#undef MTAO_TBB_ENABLED
#if defined(MTAO_TBB_ENABLED)
    tbb::parallel_for<int>(0, I.size(), [&](const int j) {
#else
    for (int j = 0; j < I.size(); ++j) {
#endif
        int& r = I(j);
        auto p = P.col(j);
        auto [grid_cell, quotient] = vertex_grid().coord(p);
        if (!g.valid_index(grid_cell)) {
#if defined(MTAO_TBB_ENABLED)
            return;
#else
            continue;
#endif
            // check if its  in an adaptive grid cell, tehn we can just use
            // that cell if the grid returns -2 we pass that through
        }
        if (r = g(grid_cell); r >= 0) {
#if defined(MTAO_TBB_ENABLED)
            return;
#else
            continue;
#endif
        }
        if (r != -1) {
            spdlog::error("Grid cell ({}) has unknown negative index {}",
                          fmt::join(grid_cell, ","), r);
#if defined(MTAO_TBB_ENABLED)
            return;
#else
            continue;
#endif
        } else if (const auto it = grid_to_cells.find(grid_cell);
                   it == grid_to_cells.end()) {
            spdlog::error("Grid cell ({}) should have cut-cells but has
none", fmt::join(grid_cell, ",")); #if defined(MTAO_TBB_ENABLED) return;
#else
            continue;
#endif
        } else {
            const auto& cell_indices = it->second;
            for (auto&& ci : cell_indices) {
                const auto& cell = cells().at(ci);
                bool contains = cell.contains(*VP, faces(), p);
                // double solid_angle = cell.solid_angle(*VP, faces(), p);
                if (contains) {
                    r = ci;
                }
            }
            if (!quiet_failures && r < 0) {
                // TODO
                // if I haven't returned yet then either this ccm is bad
                // or i'm too close to an edge. lets not assume that for now

                spdlog::warn(
                    "There are {} cells in grid cell ({}) but Mandoline "
                    "failed to find point ({}) in it. Falling back to max "
                    "solid angles",
                    cell_indices.size(), fmt::join(grid_cell, ","),
                    fmt::join(p, ","));
                int max_cell = -1;
                double max_solid_angle = std::numeric_limits<double>::min();
                for (auto&& ci : cell_indices) {
                    const auto& cell = cells().at(ci);
                    double solid_angle = cell.solid_angle(*VP, faces(), p);
                    bool contains = cell.contains(*VP, faces(), p);
                    if (solid_angle > max_solid_angle) {
                        max_cell = ci;
                        max_solid_angle = solid_angle;
                    }
                    spdlog::info(
                        "Cell not found debug: got a solid angle of {} "
                        "from "
                        "cell "
                        "{} with contains = {}",
                        solid_angle, ci, contains);
                }
                if (max_solid_angle > .5) {
                    r = max_cell;
                }
            }
        }
#if defined(MTAO_TBB_ENABLED)
    });
#else
    }
#endif
    return I;
    */
}
/*

Eigen::Ref<const ColVecs> P) const {
    const auto g = exterior_grid().cell_ownership_grid();
    const auto grid_to_cells = cut_cells_per_grid_cell();
    auto cells = get_cell_indices(P);
#if defined(MTAO_TBB_ENABLED)
    tbb::parallel_for<int>(0, cells.size(), [&](const int j) {
#else
    for (int j = 0; j < I.size(); ++j) {
#endif
        int& cell = cells(j);
        if (cell < 0) {
            cell = get_nearest_cell_index(P.col(j));
        }
#if defined(MTAO_TBB_ENABLED)
    });
#else
    }
#endif
    return cells;
}
*/

}  // namespace mandoline::operators

