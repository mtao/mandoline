

#include "mandoline/tools/cutmesh_info.hpp"

#include <mandoline/operators/boundary3.hpp>
namespace mandoline::tools {
void print_file_info(const std::string& filename, nlohmann::json* json) {
    CutCellMesh<3> ccm = CutCellMesh<3>::from_proto(filename);
    if (ccm.empty()) {
        std::cout << "Empty cutmesh or not a cutmesh!" << std::endl;
        return;
    } else {
        if (json) {
            (*json)["filename"] = filename;
        } else {
            std::cout << "Filename: " << filename << std::endl;
        }
        print_all_info(ccm, json);
    }
}

void print_all_info(const CutCellMesh<3>& ccm, nlohmann::json* json) {
    print_general_info(ccm, json);
    print_cell_info(ccm, json);
    print_face_info(ccm, json);
    print_region_info(ccm, json);
    print_volume_info(ccm, json);
    print_face_adjacent_grid_info(ccm, json);
    print_boundary_info(ccm, json);
}

void print_general_info(const CutCellMesh<3>& ccm, nlohmann::json* json) {
    if (!json) {
        std::cout << "General Cutcell Information" << std::endl;
        std::cout << "===========================" << std::endl;
    }

    {
        auto s = ccm.cell_shape();
        if (json) {
            (*json)["shape"] = {{"x", s[0]}, {"y", s[1]}, {"z", s[2]}};
        } else {
            std::cout << "Grid cell shape: ";
            std::copy(s.begin(), s.end(),
                      std::ostream_iterator<int>(std::cout, " "));
            std::cout << std::endl;
        }
    }
    if (json) {
            (*json)["vertex_count"] = ccm.origV().cols();
            (*json)["face_count"] = ccm.origF().cols();
            (*json)["cut_vertex_count"] = ccm.vertex_size();
        (*json)["cell_info"] = {
            {"total", ccm.cell_size()},
            {"cut_count", ccm.cut_cell_size()},
            {"cubic_count", ccm.exterior_grid().num_cells()}};
    } else {
        if (size_t size = ccm.cell_size(); size > 0) {
            std::cout << "Number of cells: ";
            std::cout << ccm.cell_size() << " (";
            if (size_t size = ccm.cut_cell_size(); size > 0) {
                std::cout << size << " cut-cells, ";
            }
            if (size_t size = ccm.exterior_grid().num_cells(); size > 0) {
                std::cout << size << " cubic/adaptive-cells";
            }
            std::cout << ")" << std::endl;
        } else {
            std::cout << "No cutcells!" << std::endl;
        }
    }
}

void print_face_info(const CutCellMesh<3>& ccm, nlohmann::json* json) {
    if (!json) {
        std::cout << "Face Information" << std::endl;
        std::cout << "================" << std::endl;
    }
    if (json) {
        (*json)["face_info"] = {
            {"total", ccm.face_size()},
            {"cut_count", ccm.cut_face_size()},
            {"cached_triangulations", ccm.has_triangulated_faces_cached()},
            {"cubic_count", ccm.exterior_grid().num_faces()}};
    } else {
        if (size_t size = ccm.face_size(); size > 0) {
            std::cout << "Number of faces: ";
            std::cout << ccm.face_size() << " (";

            if (size_t size = ccm.cut_face_size(); size > 0) {
                std::cout << size << " cut-faces, ";
            }
            if (size_t size = ccm.exterior_grid().num_faces(); size > 0) {
                std::cout << size << " cubic/adaptive-faces, ";
            }

            if (ccm.has_triangulated_faces_cached()) {
                std::cout << "has faces cached";
            } else {
                std::cout << "no faces cached";
            }
            std::cout << ")" << std::endl;
        } else {
            std::cout << "No cutfaces!" << std::endl;
        }
        std::cout << std::endl;
    }
}
void print_region_info(const CutCellMesh<3>& ccm, nlohmann::json* json) {
    if (!json) {
        std::cout << "Region Information" << std::endl;
        std::cout << "==================" << std::endl;
    }
    auto r = ccm.regions();
    std::map<int, int> region_counts;
    for (auto&& v : r) {
        region_counts[v]++;
    }
    if (json) {
        (*json)["cell_regions"] = r;
        (*json)["region_cell_counts"] = region_counts;
    } else {
        std::cout << region_counts.size() << " regions found" << std::endl;
        for (auto&& [id, count] : region_counts) {
            std::cout << ">>Region " << id << " has " << count << " cells"
                      << std::endl;
        }
    }
}

void print_volume_info(const CutCellMesh<3>& ccm, nlohmann::json* json) {
    auto v = ccm.cell_volumes();
    if (json) {
        (*json)["cell_volumes"] = v;
    }
}

void print_cell_info(const CutCellMesh<3>& ccm, nlohmann::json* json) {
    std::cout << "Cell Information" << std::endl;
    std::cout << "================" << std::endl;
    for (auto&& c : ccm.cells()) {
        auto g = c.grid_cell;
        for (auto&& [dim, i, m] :
             mtao::iterator::enumerate(g, ccm.cell_shape())) {
            if (i < 0 || i >= m) {
                std::cout << "Cell outside grid (dim=" << dim << ")" << i << "/"
                          << m << std::endl;
                std::cout << "Face/sign pairs: ";
                for (auto&& [fidx, b] : c) {
                    std::cout << fidx << ":" << b << " ";
                }
                std::cout << std::endl;
            }
        }
    }
}

void print_boundary_info(const CutCellMesh<3>& ccm, nlohmann::json* json) {
    std::vector<std::set<std::tuple<int, double>>> boundary(ccm.cell_size());
    auto B = operators::boundary(ccm, true);
    for (int o = 0; o < B.outerSize(); ++o) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(B, o); it; ++it) {
            boundary[it.col()].emplace(it.row(), it.value());
        }
    }

    (*json)["boundary_matrix"] = boundary;
}

void print_face_adjacent_grid_info(const CutCellMesh<3>& ccm,
                                   nlohmann::json* json) {
    std::vector<std::set<std::array<int, 3>>> face_neighbor_cells(
        ccm.face_size());
    auto B = operators::boundary(ccm, true);
    for (int o = 0; o < B.outerSize(); ++o) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(B, o); it; ++it) {
            int face_index = it.row();
            int cell_index = it.col();
            if (ccm.is_cut_cell(cell_index)) {
                const auto& c = ccm.cells()[cell_index];
                face_neighbor_cells[face_index].emplace(c.grid_cell);
            } else {
                const auto& cell = ccm.exterior_grid().cells().at(cell_index);
                auto c = cell.corner();
                auto r = cell.width();
                for (int j = c[0]; j < c[0] + r; ++j) {
                    for (int k = c[1]; k < c[1] + r; ++k) {
                        for (int l = c[2]; l < c[2] + r; ++l) {
                            face_neighbor_cells[face_index].emplace(
                                std::array<int, 3>{{j, k, l}});
                        }
                    }
                }
            }
        }
    }

    (*json)["face_adjacent_grid_cells"] = face_neighbor_cells;
}
}  // namespace mandoline::tools
