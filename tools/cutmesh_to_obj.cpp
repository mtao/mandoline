#include <fmt/format.h>

#include <mtao/cmdline_parser.hpp>
#include <mtao/logging/logger.hpp>

#include "mandoline/mesh3.hpp"
using namespace mtao::logging;
using namespace mandoline;

void write_obj(const CutCellMesh<3>& ccm, const std::string& prefix,
               const std::set<int>& indices,
               const std::string& name_suffix = ".obj") {
    std::map<int, bool> used_faces;  // bool is whether to flip
    {  // build up the faces that will be used in the obj
        auto used_cells = Eigen::VectorXd::Zero(ccm.num_cells()).eval();

        for (auto&& i : indices) {
            used_cells(i) = 1;
        }
        used_cells.setConstant(1);
        mtao::VecXd used_face_signs = ccm.boundary(true) * used_cells;
        std::cout << used_face_signs.transpose();
        std::cout << "Number of face signs available vs total face count: "
                  << used_face_signs.size() << " / " << ccm.num_faces()
                  << " = (" << ccm.num_cut_faces() << " + "
                  << ccm.exterior_grid().num_faces() << ")" << std::endl;
        // used_face_signs.setConstant(1);

        for (int i = 0; i < used_face_signs.size(); ++i) {
            if (used_face_signs(i) != 0) {
                /*
                if(ccm.is_cut_face(i)) {
                    auto&& f = ccm.cut_face(i);
                    if(f.external_boundary) {
                        auto [a,sgn]= *f.external_boundary;
                        std::cout << std::string(f) << std::endl;
                    } else {
                        std::cout << "Non external boundary boundary failure!"
                << std::endl;
                    }
                } else {
                    auto&& f = ccm.exterior_grid().face(i-ccm.num_cut_faces());
                    auto [a,b,c] = f.corner();
                    auto d = f.axis();
                    auto [u,v] = f.dual_edge;
                    std::cout << "External entry: " << a << "," << b << "," << c
                << " axis=" << d << "Dual edgE:" << u << ":" << v << std::endl;
                }
                    */
                used_faces[i] = used_face_signs(i) > 0;
            }
        }
    }
    // used_faces.clear();
    // for(int i = 0; i < ccm.exterior_grid().num_faces(); ++i) {
    //    used_faces[i+ccm.num_cut_faces()] = true;
    //}

    auto V = ccm.vertices();
    auto subV = ccm.compute_subVs();

    std::vector<mtao::ColVecs3i> Fs(used_faces.size());
    int vertex_offset = V.cols();
    for (auto&& [F, pr] : mtao::iterator::zip(Fs, used_faces)) {
        auto [face_index, flip] = pr;
        if (ccm.is_cut_face(face_index)) {
            auto&& f = ccm.cut_face(face_index);
            if (f.triangulation) {
                F = *f.triangulation;
            } else {
                F = f.triangulate(subV);
            }
            if (flip) {
                auto R = F.row(0).eval();
                F.row(0) = F.row(1);
                F.row(1) = R;
            }
        } else {
            F = ccm.exterior_grid().triangulated_face(
                face_index, ccm.num_cut_faces(), flip);
        }
    }

    std::ofstream ofs(prefix + name_suffix);

    std::tie(V, Fs) = mtao::geometry::mesh::compactify(V, Fs);
    auto F = balsa::eigen::hstack_iter(Fs.begin(), Fs.end());
    for (int i = 0; i < V.cols(); ++i) {
        auto v = V.col(i);
        ofs << "v " << v.transpose() << std::endl;
    }

    for (int i = 0; i < F.cols(); ++i) {
        auto f = F.col(i).array() + 1;
        ofs << "f " << f.transpose() << std::endl;
    }
}
void write_obj(const CutCellMesh<3>& ccm, const std::string& filename) {
    std::set<int> inds;
    for (int i = 0; i < ccm.num_cells(); ++i) {
        inds.insert(i);
    }
    write_obj(ccm, filename, inds);
}

void write_obj_regions(const CutCellMesh<3>& ccm, const std::string& filename) {
    std::vector<int> R = ccm.regions();
    std::map<int, std::set<int>> Rs;
    for (auto&& [i, r] : mtao::iterator::enumerate(R)) {
        Rs[r].insert(i);
    }
    for (auto&& [i, inds] : Rs) {
        write_obj(ccm, filename, inds, fmt::format("-r{}.obj", i));
    }
}
void write_mesh_obj_separate(const CutCellMesh<3>& ccm,
                             const std::string& prefix) {
    auto R = ccm.regions();
    int i = 0;
#pragma omp parallel for
    for (i = 0; i < ccm.num_cells(); ++i) {
        write_obj(ccm, prefix, {i}, fmt::format("-{}-r{}.obj", i, R[i]));
    }
}
void write_obj_separate(const CutCellMesh<3>& ccm, const std::string& prefix) {
    auto R = ccm.regions();
    int i = 0;
#pragma omp parallel for
    for (i = 0; i < ccm.num_cells(); ++i) {
        // write_obj(prefix,i);
        write_obj(ccm, prefix, {i}, fmt::format("-{}.obj", i));
    }
}

int main(int argc, char* argv[]) {
    auto&& log = make_logger("profiler", mtao::logging::Level::All);
    mtao::CommandLineParser clp;
    clp.add_option("open-regions", false);
    clp.add_option("write-regions", true);
    clp.add_option("write-monolithic", true);
    clp.add_option("write-separate", true);
    clp.add_option("normalize", false);
    clp.add_option("normalize-unit", false);
    clp.add_option("cell-grid-ownership", false);
    clp.parse(argc, argv);

    if (clp.args().size() < 1) {
        fatal() << "No input mesh filename!";
        return {};
    }
    if (clp.args().size() < 2) {
        fatal() << "No output mesh filename!";
    }

    std::string input_cutmesh = clp.arg(0);
    std::string output_prefix = clp.arg(1);
    bool write_regions = clp.optT<bool>("write-regions");
    bool normalize = clp.optT<bool>("normalize");
    bool normalize_gc = clp.optT<bool>("normalize-unit");
    bool write_monolithic = clp.optT<bool>("write-monolithic");
    bool write_separate = clp.optT<bool>("write-separate");
    bool open_regions = clp.optT<bool>("open-regions");
    bool cell_grid_ownership = clp.optT<bool>("cell-grid-ownership");

    CutCellMesh<3> ccm = CutCellMesh<3>::from_proto(input_cutmesh);

    ccm.triangulate_faces(false);

    std::set<int> inds;
    for (int i = 0; i < ccm.num_cells(); ++i) {
        if (!ccm.is_cut_cell(i)) {
            inds.insert(i);
        }
    }
    write_obj(ccm, "exterior", inds);
    return 0;

    if (write_regions) {
        write_obj_regions(ccm, output_prefix);
    }
    if (write_monolithic) {
        write_obj(ccm, output_prefix);
    }
    if (write_separate) {
        write_obj_separate(ccm, output_prefix);
    }
    if (cell_grid_ownership) {
        std::ofstream ofs(output_prefix + "_cells.txt");
        auto C = ccm.cell_centroids();
        for (auto&& [i, cell] : mtao::iterator::enumerate(ccm.cells())) {
            ofs << i << " ";
            auto cent = C.col(i);
            auto [c, q] = ccm.coord(cent);
            std::copy(c.begin(), c.end(), std::ostream_iterator<int>(ofs, " "));
            ofs << std::endl;
        }
        for (auto&& [i, c] : ccm.exterior_grid().cells()) {
            auto center = c.corner();
            int w = c.width();
            for (auto&& c : center) {
                c += w / 2;
            }
            ofs << i << " ";
            std::copy(center.begin(), center.end(),
                      std::ostream_iterator<int>(ofs, " "));
            ofs << std::endl;
        }
    }
    return 0;
}

