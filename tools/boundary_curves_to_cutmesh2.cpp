#include <algorithm>
#include <iterator>
#include <mtao/types.hpp>

#include <mandoline/construction/generator2.hpp>
#include <cxxopts.hpp>
#include <spdlog/spdlog.h>


// Input:
// Vertices
// vertex_idx, xcoord, ycoord
//
// Edges
// edge_idx, vertex_idx1, vertex_idx2, left_element_idx, right_element_idx
std::tuple<mtao::ColVecs2d, mtao::ColVecs4i> read_curve_file(const std::string &filename) {
    std::ifstream ifs(filename);
    std::string line;

    enum class State { Init,
                       Vertices,
                       Edges };
    State state = State::Init;

    auto trim = [](std::string &str) {
        str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](std::string::value_type v) {
                      return !std::isspace(v);
                  }));
        str.erase(std::find_if(str.rbegin(),str.rend(),[](std::string::value_type v) {
                return !std::isspace(v);
                }).base(), str.end());
    };

    auto tolower = [](std::string &str) {
        std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::tolower(c); });
    };

    using VertexEntry = std::array<double, 2>;
    using EdgeEntry = std::array<int, 4>;

    std::map<int, VertexEntry> vertices;
    std::map<int, EdgeEntry> edges;

    auto process_vertex = [&](const std::vector<std::string>& toks) {
        if (toks.size() < std::tuple_size_v<VertexEntry>) {
            spdlog::error("Vertex with wrong number of tokens! ignoring");
            return;
        }
        VertexEntry e;

        int idx = std::stoi(toks[0]);
        e[0] = std::stod(toks[1]);
        e[1] = std::stod(toks[2]);
        vertices[idx] = e;

    };

    auto process_edge = [&](const std::vector<std::string> &toks) {
        if (toks.size() < std::tuple_size_v<EdgeEntry>) {
            if (toks.size() != 3) {
                spdlog::error("Edge with wrong number of tokens! ignoring");
                return;
            }
        }

        EdgeEntry e;
        int idx = std::stoi(toks[0]);
        e[0] = std::stoi(toks[1]);
        e[1] = std::stoi(toks[2]);
        if (toks.size() == 3) {
            e[2] = -1;
            e[3] = -1;
        } else {
            e[2] = std::stoi(toks[3]);
            e[3] = std::stoi(toks[4]);
        }
        edges[idx] = e;
    };

    while (getline(ifs, line)) {
        if (line.empty()) continue;
        // be lazy about csv files and turn them into space separated values
        std::transform(line.begin(), line.end(), line.begin(), [](unsigned char c) -> unsigned char {
            if (c == ',') { return ' '; }
            return c;
        });

        // split into tokens
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_insert_iterator<std::vector<std::string>>(tokens));

        // if this was all whitespace flff move on
        if (tokens.empty()) { continue; }

        // maybe its something like "vertices" or "edges"
        if (tokens.size() == 1) {
            trim(line);
            tolower(line);
            if (line == "vertices") {
                state = State::Vertices;
                continue;
            } else if (line == "edges") {
                state = State::Edges;
            } else if (line == "elements") {// if we're starting to read elements this is something that was already cut. lets just cut it again!
                break;
            }
        } else {
            try {
                switch (state) {
                case State::Vertices:
                    process_vertex(tokens);
                case State::Edges:
                    process_edge(tokens);
                case State::Init:
                    continue;
                }
            } catch (const std::invalid_argument &e) {
                continue;
            } catch (const std::out_of_range &e) {
                continue;
            }
        }
    }


    int max_vidx = std::max_element(vertices.begin(), vertices.end())->first;
    int max_eidx = std::max_element(edges.begin(), edges.end())->first;
    mtao::ColVecs2d V(2, max_vidx + 1);
    mtao::ColVecs4i E(4, max_vidx + 1);
    V.setZero();
    E.setZero();
    for (auto [k, v] : vertices) {
        V.col(k) = mtao::eigen::stl2eigen(v);
    }
    for (auto [k, v] : edges) {
        E.col(k) = mtao::eigen::stl2eigen(v);
    }
    return { V, E };
}

//pass in the inputE so we can extract this left/right element idx stuff
void write_boundary_curve(const mandoline::CutCellMesh<2> &ccm, const mtao::ColVecs4i &inputE, const std::string &filename) {

    std::ofstream ofs(filename);
    ofs << "Vertices" << std::endl;
    auto V = ccm.vertices();
    for (int i = 0; i < V.cols(); ++i) {
        auto v = V.col(i);
        ofs << i << ", " << v(0) << ", " << v(1) << "\n";
    }

    ofs << std::endl;
    ofs << "Edges" << std::endl;

    for (auto &&[idx, ce] : mtao::iterator::enumerate(ccm.cut_edges())) {
        ofs << idx << ", " << ce.indices[0] << ", " << ce.indices[1] << ", ";
        if (ce.is_mesh_edge()) {
            auto &&ie = ccm.mesh_edges().at(idx);
            auto pe = inputE.col(ce.as_edge_id());
            if (ie.ts(0) > ie.ts(1)) {
                ofs << pe(3) << ", " << pe(2) << "\n";
            } else {
                ofs << pe(2) << ", " << pe(3) << "\n";
            }
        } else {
            ofs << -1 << " " << -1 << "\n";
        }
    }

    {
        int offset = ccm.cut_edges().size();
    auto&& eg = ccm.exterior_grid;
        for (auto &&[idx, ce] : mtao::iterator::enumerate(eg.boundary_facet_pairs())) {
            
            int axis = eg.get_face_axis(idx);
            std::array<int,2> e;
            std::array<int,2> corner;
            if(ce[1] >= 0) {
                corner = eg.cell_coord(ce[1]);

            } else if(ce[0] >= 0) {
                corner = eg.cell_coord(ce[0]);
                corner[axis]+=1;
            } else {
                spdlog::error("File a bug report! mesh has an exterior grid entry with nothing in it!");
            }
            e[0] = eg.Base::vertex_index(corner);
            corner[1-axis]+=1;
            e[1] = eg.Base::vertex_index(corner);

            ofs << (idx + offset) << ", " << e[0] << ", " << e[1] << ", ";
            ofs << -1 << " " << -1 << "\n";
        }
    }
    ofs << std::endl;
    ofs << "Elements" << std::endl;

    for (auto &&[idx, face] : mtao::iterator::enumerate(ccm.cut_faces())) {
        for(auto&& indices: face.indices) {
            ofs << idx << ", " << indices.size() << ", ";
            std::copy(indices.begin(), indices.end(), std::ostream_iterator<int>(ofs, ", "));
            ofs << " " << face.region;
            ofs << "\n";
        }
    }
    {
    int offset = ccm.cut_faces().size();
    auto&& eg = ccm.exterior_grid;
    for(int i = 0; i < eg.num_cells(); ++i) {
        auto cc = eg.cell_coord(i);
        auto&& r = eg.region(i);

        ofs << (offset + i) << ", " << 4 <<", ";
        ofs << eg.Base::vertex_index(cc) << ", ";
        cc[1]++;
        ofs << eg.Base::vertex_index(cc) << ", ";
        cc[0]++;
        ofs << eg.Base::vertex_index(cc) << ", ";
        cc[1]--;
        ofs << eg.Base::vertex_index(cc) << ", ";
        ofs << r << std::endl;
    }
    ofs << std::flush;
    }
}


int main(int argc, char *argv[]) {
    cxxopts::Options options("boundary_curve_to_cutmesh2", "A brief description");

    options.add_options()("input", "input", cxxopts::value<std::string>())(
      "output", "output mesh file", cxxopts::value<std::string>())(
      "Nx", "Grid dimension in x", cxxopts::value<int>()->default_value("5"))(
      "Ny", "Grid dimension in y", cxxopts::value<int>()->default_value("5"))
      ("mx", "min value in x", cxxopts::value<double>()->default_value("0"))("my", "min value in y", cxxopts::value<double>()->default_value("0"))("Mx", "max value in x", cxxopts::value<double>()->default_value("10"))("My", "max value in y", cxxopts::value<double>()->default_value("10"))(
      "h,help", "Print usage");
    options.parse_positional({ "input", "output" });
    options.positional_help({ "<input file> <output_file>" });

    auto result = options.parse(argc, argv);
    Eigen::AlignedBox<double, 2> bbox;
    bbox.min()
      << result["mx"].as<double>()
      , result["my"].as<double>();
    bbox.max()
      << result["Mx"].as<double>()
      , result["My"].as<double>();

    int Nx = result["Nx"].as<int>();
    int Ny = result["Ny"].as<int>();
    auto grid = mtao::geometry::grid::StaggeredGrid2d::from_bbox(bbox, { { Nx, Ny } });
    auto [oV, oE] = read_curve_file(result["input"].as<std::string>());




    mandoline::construction::CutCellGenerator<2> ccg(oV, grid);
    ccg.add_boundary_elements(oE.topRows<2>());
    ccg.bake();

    auto ccm = ccg.generate();
    write_boundary_curve(ccm, oE, result["output"].as<std::string>());
}

// OUTPUT:
// Vertices
// vertex_idx, xcoord, ycoord
//
// Edges
// edge_idx, vertex_idx1, vertex_idx2, left_element_idx, right_element_idx
//
// Elements
// element_idx, num_vertices, vertex_idx1, vertex_idx2, ..., vertex_idxn, componentID
