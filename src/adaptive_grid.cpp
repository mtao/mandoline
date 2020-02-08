#include "mandoline/adaptive_grid.hpp"
#include <mtao/eigen/stl2eigen.hpp>
#include <set>
#include "mandoline/proto_util.hpp"
#include <iterator>
namespace mandoline {
    template <typename GridB>
        void print_gridb(const GridB& g) {

            if constexpr(GridB::D == 3) 
            {
                for(int i = 0; i < g.shape(0); ++i) {
                    for(int j = 0; j < g.shape(1); ++j) {
                        for(int k = 0; k < g.shape(2); ++k) {
                            if(g(i,j,k)) {
                                std::cout << "o";
                            } else {
                                std::cout << ".";
                            }
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }

        }
    template <typename GridB>
        void print_grid(const GridB& g) {

            if constexpr(GridB::D == 3) 
            {
                for(int i = 0; i < g.shape(0); ++i) {
                    for(int j = 0; j < g.shape(1); ++j) {
                        for(int k = 0; k < g.shape(2); ++k) {
                            std::cout << g(i,j,k) << " ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }

        }

    void   AdaptiveGrid::Cell::serialize(protobuf::Cube& c) const {
        protobuf::serialize(corner(),*c.mutable_corner());
        c.set_radius(width());

    }
    AdaptiveGrid::Cell AdaptiveGrid::Cell::from_proto(const protobuf::Cube& c) {
        Cell ret;
        protobuf::deserialize(c.corner(),std::get<0>(ret));
        std::get<1>(ret) = c.radius();
        return ret;

    }
    bool AdaptiveGrid::Cell::is_inside(const Vec& p) const {
        Eigen::Map<const mtao::Vec3i> c(corner().data());
        return (c.array().cast<double>() <= p.array()).all() 
            && (p.array() < (c.array()+width()).cast<double>()).all();
    }
    void   AdaptiveGrid::Square::serialize(protobuf::Square& c) const {
        protobuf::serialize(corner(),*c.mutable_corner());
        c.set_radius(width());
        c.set_axis(axis());

    }
    AdaptiveGrid::Square AdaptiveGrid::Square::from_proto(const protobuf::Square& c) {
        Square ret;
        protobuf::deserialize(c.corner(),std::get<0>(ret));
        std::get<1>(ret) = c.radius();
        std::get<2>(ret) = c.axis();
        return ret;

    }
    void   AdaptiveGrid::Face::serialize(protobuf::Face& c) const {
        Square::serialize(*c.mutable_geometry());
        protobuf::serialize(dual_edge, *c.mutable_dual_edge());

    }
    AdaptiveGrid::Face AdaptiveGrid::Face::from_proto(const protobuf::Face& c) {
        Face ret;
        ret.Square::operator=(Square::from_proto(c.geometry()));
        protobuf::deserialize(c.dual_edge(), ret.dual_edge);
        return ret;

    }

    auto AdaptiveGrid::cell_ownership_grid() const -> GridData3i {
        return grid_from_cells(cell_shape(),m_cells);
    }
    auto AdaptiveGrid::grid_from_cells(const coord_type& shape, const std::map<int,Cell>& cells) -> GridData3i {
        GridData3i grid = GridData3i::Constant(-1,shape);
        for(auto&& [cid,desc]: cells) {
            auto&& abc = desc.corner();
            auto&& jump = desc.width();
            for(int aa = abc[0]; aa < abc[0] + jump; ++aa) {
                for(int bb = abc[1]; bb < abc[1] + jump; ++bb) {
                    for(int cc = abc[2]; cc < abc[2] + jump; ++cc) {
                        auto& c = grid(aa,bb,cc);
                        if(c == -1) {
                            c = cid;
                        } else {
                            std::cout << "Overlapping cells! Checked" << aa << "," << bb << "," << cc << " with cid " <<cid << " got " << c << "instead" << std::endl;
                        }

                    }
                }
            }
        }
        return grid;

    }
    std::array<int,4> AdaptiveGrid::face(const Cell& c, int axis, bool sign) const {
        coord_type p = c.corner();
        if(sign) {
            p[axis]++;
        }
        std::array<int,4> ret;
        ret[0] = vertex_index(p);
        int u = (axis+1)%3;
        int v = (axis+2)%3;
        p[u]++;
        ret[0] = vertex_index(p);
        p[v]++;
        ret[0] = vertex_index(p);
        p[u]--;
        ret[0] = vertex_index(p);
        return ret;
    }
    std::array<int,4> AdaptiveGrid::face(int idx, int axis, bool sign) const {
        return face(cell(idx),axis,sign);
    }
    auto AdaptiveGrid::Cell::vertex(int a, int b, int c) const -> coord_type {
        return coord_type{{
            corner()[0] + a*width(), 
                corner()[1] + b*width(), 
                corner()[2] + c*width()}};

    }
    auto AdaptiveGrid::Square::vertex(int a, int b) const -> coord_type {
        coord_type c = corner();
        int d = dimension();
        c[(d+1)%3] += width() * a;
        c[(d+2)%3] += width() * b;

        return c;
    }
    mtao::ColVecs3i AdaptiveGrid::triangulated(int idx) const {
        return triangulated(cell(idx));
    }
    mtao::ColVecs3i AdaptiveGrid::triangulated(const Cell& c) const {
        auto f = []() {
            Eigen::MatrixXi F(12,3);
            F << 0, 2, 1, //face front
              0, 3, 2,
              2, 3, 4, //face top
              2, 4, 5,
              1, 2, 5, //face right
              1, 5, 6,
              0, 7, 4, //face left
              0, 4, 3,
              5, 4, 7, //face back
              5, 7, 6,
              0, 6, 7, //face bottom
              0, 1, 6;
            return F.transpose().eval();
        };


        std::array<int,8> v = {{
            static_cast<int>(vertex_index(c.vertex(0, 0, 0))),
                static_cast<int>(vertex_index(c.vertex(1, 0, 0))),
                static_cast<int>(vertex_index(c.vertex(1, 1, 0))),
                static_cast<int>(vertex_index(c.vertex(0, 1, 0))),
                static_cast<int>(vertex_index(c.vertex(0, 1, 1))),
                static_cast<int>(vertex_index(c.vertex(1, 1, 1))),
                static_cast<int>(vertex_index(c.vertex(1, 0, 1))),
                static_cast<int>(vertex_index(c.vertex(0, 0, 1)))}};
        const static mtao::ColVecs3i F = f();
        return F.unaryExpr([&](int idx) -> int {
                return v[idx];
                });
    }

    auto AdaptiveGrid::boundary_triplets(int offset) const -> std::vector<Eigen::Triplet<double>> {
        std::vector<Eigen::Triplet<double>> trips;
        for(auto&& [i,face]: mtao::iterator::enumerate(m_faces)) {
            auto& e = face.dual_edge;
            if(is_valid_edge(e)) {
                auto [l,h] = e;
                trips.emplace_back(offset+i,l,-1);
                trips.emplace_back(offset+i,h,1);
            }
        }
        return trips;
    }
    auto AdaptiveGrid::faces(const GridData3i& grid) const -> std::vector<Face> {
        auto myless = [](const Face& a, const Face& b) {
            return std::less<Edge>()(a.dual_edge,b.dual_edge);
        };
        std::set<Face,decltype(myless)> faces(myless);
        auto add_bdry = [&](int d, int a, int b) {
            Edge dual_edge{{a,b}};
            //if(is_valid_edge(dual_edge)) {
            {
                int axis=d, width;
                coord_type corner;
                if(a == -1) {
                    auto&& c = cell(b);
                    corner = c.corner();
                    width = c.width();
                } else if(b == -1) {
                    auto&& c = cell(a);
                    corner = c.corner();
                    width = c.width();
                    corner[d] += width;
                } else {
                    auto&& ca = cell(a);
                    auto&& cb = cell(b);
                    width = std::min(ca.width(),cb.width());
                    //if higher one is the smaller one we just use it
                    if(width == cb.width()) {//checking cb is important
                        corner = cb.corner();
                    } else {
                        corner = ca.corner();
                        corner[d] += width;
                    }
                }
                faces.emplace(Square{corner,axis,width},dual_edge);
            }

        };
        for(int d = 0; d < 3; ++d) {
            auto shape = cell_shape();
            shape[d]++;
            coord_type abc;
            for(int& a = abc[0] = 0; a < shape[0]; ++a) {
                for(int& b = abc[1] = 0; b < shape[1]; ++b) {
                    for(int& c = abc[2] = 0; c < shape[2]; ++c) {
                        bool low =  abc[d] == 0;
                        bool high = abc[d] == shape[d]-1;
                        if( low ^ high) {
                            if(low) {
                                add_bdry(d,-1,grid(abc));
                            } else {
                                // if we're high we have to go down one grid cell
                                coord_type tmp = abc;
                                tmp[d]--;
                                add_bdry(d,grid(tmp),-1);
                            }
                        } else {
                            int pidx = cell_index(abc);
                            coord_type tmp = abc;
                            tmp[d]--;
                            int nidx = cell_index(tmp);
                            int pid = grid.get(pidx);
                            int nid = grid.get(nidx);
                            if(pid >= 0 && nid >= 0) {
                                add_bdry(d,nid,pid);
                            }
                        }
                    }
                }
            }
        }

        std::vector<Face> faces_vec(faces.size());
        std::copy(faces.begin(),faces.end(),faces_vec.begin());

        return faces_vec;
    }
    void AdaptiveGrid::make_faces() {
        m_faces = faces(cell_ownership_grid());
    }

    mtao::VecXd AdaptiveGrid::dual_edge_lengths() const {

        auto&dx = Base::dx();
        mtao::VecXd ret(num_faces());
        for(auto&& [i,f]: mtao::iterator::enumerate(m_faces)) {
            auto&& e = f.dual_edge;
            if(!is_valid_edge(e)) continue;
            auto [a,b] = e;
            ret(i) = (dx.asDiagonal() * (cell(a).center() - cell(b).center())).norm();
        }
        return ret;
    }
    void AdaptiveGrid::cell_centroids(mtao::ColVecs3d& R) const {
        for(auto&& [i,c]: cells()) {
            auto cent = R.col(i);
            cent = mtao::eigen::stl2eigen(c.corner()).cast<double>();
            double w = c.width();
            cent.array() += w / 2.0;
            cent = vertex_grid().world_coord(cent);
        }
    }
    mtao::VecXd AdaptiveGrid::cell_volumes() const {
        mtao::VecXd R(num_cells());
        if(num_cells() == 0) {
            return {};
        } else {
            int offset = m_cells.begin()->first;
            double vol = dx().prod();
            for(auto&& [i,c]: cells()) {
                double w = c.width();
                R(i-offset) = vol * (w*w*w);
            }
            return R;
        }
    }

    mtao::VecXd AdaptiveGrid::face_volumes() const {
        auto&dx = Base::dx();
        mtao::VecXd ret(num_faces());
        mtao::Vec3d dws = mtao::Vec3d::Ones();
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 2; ++j) {
                dws(i) = dx((j+i)%3);
            }
        }
        for(auto&& [i,f]: mtao::iterator::enumerate(m_faces)) {
            auto&& e = f.dual_edge;
            if(!is_valid_edge(e)) continue;
            int w= f.width();
            ret(i) = w * w * dws(f.axis());
        }
        return ret;
    }
    int AdaptiveGrid::num_faces() const {
        return m_faces.size();
    }
    int AdaptiveGrid::num_cells() const {
        return cells().size();
    }

    mtao::ColVecs3d AdaptiveGrid::boundary_centroids() const {
        mtao::ColVecs3d C(3,num_faces());
        for(auto&& [index,f]: mtao::iterator::enumerate(m_faces)) {
            auto&& e = f.dual_edge;
            auto cent = C.col(index);
            cent.setZero();
            if(!is_valid_edge(e)) continue;
            auto [a,b] = e;
            int k = 0;

            auto& ca = cell(a);
            auto& cb = cell(b);
            mtao::Vec3d cd = ca.center() - cb.center();
            int minwidth = std::min(ca.width(),cb.width());
            cd.cwiseAbs().maxCoeff(&k);
            coord_type corner = ca.width() < cb.width() ? ca.corner():cb.corner();

            if(cd(k) > 0) {
                corner[k] += ca.width();
            }
            {
                const double w = minwidth;
                const coord_type c = corner;
                const int u = (k+1)%3;
                const int v = (k+2)%3;
                cent = mtao::eigen::stl2eigen(c).cast<double>();
                cent(u) +=w/2.0;
                cent(v) +=w/2.0;
                cent = vertex_grid().world_coord(cent);

            }
        }
        return C;
    }

    std::vector<Eigen::Triplet<double>> AdaptiveGrid::grid_face_projection(int offset) const {
        std::vector<Eigen::Triplet<double>> trips;
        mtao::VecXd ret(num_faces());
        trips.reserve(form_size<2>());
        for(auto&& [index,f]: mtao::iterator::enumerate(m_faces)) {
            auto&& e = f.dual_edge;
            if(!is_valid_edge(e)) continue;
            auto [a,b] = e;
            int k = 0;

            auto& ca = cell(a);
            auto& cb = cell(b);
            mtao::Vec3d cd = ca.center() - cb.center();
            int minwidth = std::min(ca.width(),cb.width());
            cd.cwiseAbs().maxCoeff(&k);
            coord_type corner = ca.width() < cb.width() ? ca.corner():cb.corner();

            if(cd(k) > 0) {
                corner[k] += ca.width();
            }
            {
                const double w = minwidth;
                const coord_type c = corner;
                const int u = (k+1)%3;
                const int v = (k+2)%3;
                double vol = 1/(w*w);
                coord_type a = c;
                for(int& i = a[u] = c[u]; i < w+c[u]; ++i) {
                    for(int& j = a[v] = c[v]; j < w+c[v]; ++j) {
                        const int row = offset+index;
                        const int col = staggered_index<2>(a,k);
                        trips.emplace_back(row, col, vol);
                    }
                }
            }
        }
        return trips;
    }
    std::vector<Eigen::Triplet<double>> AdaptiveGrid::grid_cell_projection() const {
        std::vector<Eigen::Triplet<double>> trips;
        for(auto&& [col,cell]: m_cells) {
            auto&& c= cell.corner();
            auto&& w = cell.width();
            double vol = 1. / (w*w*w);

            coord_type a;
            for(int& i = a[0]= c[0]; i < w+c[0]; ++i) {
                for(int& j = a[1]= c[1]; j < w+c[1]; ++j) {
                    for(int& k  = a[2]= c[2]; k < w+c[2]; ++k) {
                        trips.emplace_back(cell_index(a), col, vol);
                    }
                }
            }
        }
        return trips;
    }

    auto AdaptiveGrid::boundary_pairs() const -> std::vector<Edge> {
        std::vector<Edge> R(m_faces.size());
        std::transform(m_faces.begin(),m_faces.end(),R.begin(),[](const Face& f) {
                return f.dual_edge;
                });
        return R;
    }
    int AdaptiveGrid::get_cell_index(const Vec& p) const {
        for(auto&& [i,c]: cells()) {
            if(c.is_inside(p)) {
                return i;
            }
        }
        return -1;
    }
    int AdaptiveGrid::num_edges() const {
        return m_edges.size();
    }
    mtao::ColVecs2i AdaptiveGrid::edges() const {
        return mtao::eigen::stl2eigen(m_edges); 
    }
}
