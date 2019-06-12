#include "mandoline/adaptive_grid.hpp"
#include <mtao/eigen/stl2eigen.hpp>
#include <set>
#include "mandoline/proto_util.hpp"
#include <iterator>
namespace mandoline {

    void   AdaptiveGrid::Cell::serialize(CutMeshProto::Cube& c) const {
        protobuf::serialize(corner(),*c.mutable_corner());
        c.set_radius(width());

    }
    AdaptiveGrid::Cell AdaptiveGrid::Cell::from_proto(const CutMeshProto::Cube& c) {
        Cell ret;
        protobuf::deserialize(c.corner(),std::get<0>(ret));
        std::get<1>(ret) = c.radius();
        return ret;

    }

    const AdaptiveGridFactory::Indexer AdaptiveGridFactory::cmask_indexer(std::array<int,3>{{width,width,width}});
    const AdaptiveGridFactory::Indexer AdaptiveGridFactory::vmask_indexer(std::array<int,3>{{1+width,1+width,1+width}});


    std::array<AdaptiveGridFactory::Indexer,3> make_edge_indexers() {
        constexpr static int width = AdaptiveGridFactory::width;

        std::array<AdaptiveGridFactory::Indexer,3> ret;
        for(auto&& [k,g]: mtao::iterator::enumerate(ret)){ 
            std::array<int,3> shape{{width,width,width}};
            for(auto&& [i,s]: mtao::iterator::enumerate(shape)) { if(i != k) s++; }
            g = AdaptiveGridFactory::Indexer(shape);
        }
        return ret;
    }
    std::array<AdaptiveGridFactory::coord_type,3> make_mask_edge_shapes() {
        constexpr static int width = AdaptiveGridFactory::width;

        std::array<AdaptiveGridFactory::coord_type,3> ret;
        for(auto&& [k,g]: mtao::iterator::enumerate(ret)){ 
            std::array<int,3> shape{{width,width,width}};
            for(auto&& [i,s]: mtao::iterator::enumerate(shape)) { if(i != k) s++; }
            g = shape;
        }
        return ret;
    }
    const std::array<AdaptiveGridFactory::Indexer,3> AdaptiveGridFactory::mask_edge_indexers = make_edge_indexers();
    const std::array<AdaptiveGridFactory::coord_type,3> AdaptiveGridFactory::mask_edge_shapes = make_mask_edge_shapes();

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
    void AdaptiveGridFactory::make_edges(const std::optional<int>& max_level) {
        auto es = compute_edges(max_level);
        edges = mtao::eigen::stl2eigen(std::get<0>(es));
        //boundary_edges = mtao::eigen::stl2eigen(es[1]);
    }
    auto AdaptiveGridFactory::compute_axial_edges(const std::optional<int>& max_level) const -> std::tuple<std::array<std::set<Edge>,3>,AxialBEdgeMap> {

        std::cout << "Computing axial edges" << std::endl;
        std::array<coord_type,3> edge_shapes{{original.shape(),original.shape(),original.shape()}};
        edge_shapes[0][1]++;
        edge_shapes[0][2]++;
        edge_shapes[1][0]++;
        edge_shapes[1][2]++;
        edge_shapes[2][0]++;
        edge_shapes[2][1]++;
        std::array<GridData3i,3> edge_coverage = {{
            GridData3i::Constant(0,edge_shapes[0]),
                GridData3i::Constant(0,edge_shapes[1]),
                GridData3i::Constant(0,edge_shapes[2])
        }};
        std::array<std::set<Edge>,3> edges;
        AxialBEdgeMap bedges;
        for(auto&& [level,g]: mtao::iterator::enumerate(levels)) {
            /*
               std::array<GridData3i,3> edge_coverage = {{
               GridData3i::Constant(0,edge_shapes[0]),
               GridData3i::Constant(0,edge_shapes[1]),
               GridData3i::Constant(0,edge_shapes[2])
               }};
               for(auto&& [dim,E]: mtao::iterator::enumerate(edges)) {
               for(auto&& [ai,bi]: E) {
               auto a = vertex_indexer.unindex(ai);
               auto b = vertex_indexer.unindex(bi);
               for(auto abc = a; abc[dim] < b[dim]; ++abc[dim]) {
               edge_coverage[dim](abc)++;
               }
               }
               }
               for(auto&& [a,b]: mtao::iterator::enumerate(edge_coverage)) {
               std::cout << "Edge coverage " << a << ": " << std::endl;
               print_grid(b);
               }
               */

            if(max_level && *max_level == level) {
                if(level == 0) {

                    auto [Es,BEs] = make_edges(original, 0);
                    for(auto&& [A,B,C,D]: mtao::iterator::zip(edges,Es,bedges,BEs)) {
                        A.insert(B.begin(),B.end());
                        C.insert(D.begin(),D.end());
                    }
                    break;
                } else {
                    auto& g = levels[level-1];
                    GridData3 A(g.shape());
                    std::transform(g.begin(),g.end(),A.begin(),[&](const ActiveMask& mask) -> bool {
                            return mask.any();
                            });
                    auto [Es,BEs] = make_edges(A, level);
                    for(auto&& [A,B,C,D]: mtao::iterator::zip(edges,Es,bedges,BEs)) {
                        A.insert(B.begin(),B.end());
                        C.insert(D.begin(),D.end());
                    }
                    break;
                }
            } else {
                auto [Es,BEs] = make_edges(g, level);
                for(auto&& [A,B,C,D]: mtao::iterator::zip(edges,Es,bedges,BEs)) {
                    A.insert(B.begin(),B.end());
                    C.insert(D.begin(),D.end());
                }


            }

        }
        /*
           for(auto&& [dim,E]: mtao::iterator::enumerate(edges)) {
           for(auto&& [ai,bi]: E) {
           auto a = vertex_indexer.unindex(ai);
           auto b = vertex_indexer.unindex(bi);
           for(auto abc = a; abc[dim] < b[dim]; ++abc[dim]) {
           edge_coverage[dim](abc)++;
           }
           }
           }
           for(auto&& [a,b]: mtao::iterator::enumerate(edge_coverage)) {
           std::cout << "Edge coverage " << a << ": " << std::endl;
           print_grid(b);
           }
           */

        return {edges,bedges};

    }
    template <typename GridAccessor>
        void AdaptiveGridFactory::set_edge_masks(const coord_type& shape, const coord_type& corner, std::array<GridData3,3>& edge_masks, const GridAccessor& accessor) const {
            for(int a = 0; a < shape[0]; ++a) {
                for(int b = 0; b < shape[1]; ++b) {
                    for(int c = 0; c < shape[2]; ++c) {
                        coord_type abc{{a,b,c}};

                        if(accessor(abc)) {

                            mtao::eigen::stl2eigen(abc) += mtao::eigen::stl2eigen(corner);
                            for(int dim = 0; dim < 3; ++dim) {
                                int x = (dim+1)%3;
                                int y = (dim+2)%3;
                                for(int i = 0; i < 2; ++i) {
                                    for(int j = 0; j < 2; ++j) {
                                        coord_type c = abc;
                                        c[x]+=i;
                                        c[y]+=j;
                                        edge_masks[dim](c) = false;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    auto AdaptiveGridFactory::compute_edges(const std::optional<int>& max_level) const -> std::tuple<std::set<Edge>,AxialBEdgeMap> {
        std::set<Edge> edges;
        auto [Es,BEs] = compute_axial_edges(max_level);
        for(auto&& E: Es) {
            edges.insert(E.begin(),E.end());
        }
        return {edges,BEs};
    }
    auto AdaptiveGridFactory::get_edges(const std::array<GridData3,3>& edge_masks, int level, const coord_type& offset) const -> std::tuple<std::array<std::set<Edge>,3>,AxialBEdgeMap> {
        int jump = get_jump(level);
        const auto offset_map = mtao::eigen::stl2eigen(offset);
        std::array<std::set<Edge>,3> edges;
        AxialBEdgeMap bedges;
        std::bitset<3> boundary;
        auto& mask = levels_mask[level];
        auto eorigmasks = make_edge_mask(mask,level);
        for(auto&& [dim,g,og]: mtao::iterator::enumerate(edge_masks,eorigmasks)) {
            int bidx0 = (dim+1)%3;
            int bidx1 = (dim+2)%3;
            for(int a = 0; a < g.shape()[0]; ++a) {
                boundary[0] = (a == 0 || a == g.shape()[0]-1);
                for(int b = 0; b < g.shape()[1]; ++b) {
                    boundary[1] = (b == 0 || b == g.shape()[1]-1);
                    for(int c = 0; c < g.shape()[2]; ++c) {
                        boundary[2] = (c == 0 || c == g.shape()[2]-1);
                        coord_type abc{{a,b,c}};
                        if(g(abc)) {
                            coord_type abc2;
                            mtao::eigen::stl2eigen(abc2) = jump * mtao::eigen::stl2eigen(abc) + offset_map;
                            auto e = get_edge(abc2,jump,dim);
                            edges[dim].emplace(e);
                            bool is_boundary = boundary[bidx0] || boundary[bidx1];
                            if(is_boundary) {
                                if(boundary[bidx0]) {
                                    bedges[bidx0][abc[bidx0]].insert(e);
                                }
                                if(boundary[bidx0]) {
                                    bedges[bidx1][abc[bidx1]].insert(e);
                                }
                            } else {
                                for(int j = 0; j < 2 && !is_boundary; ++j) {
                                    int d = (dim+j+1)%3;
                                    for(int k = 0; k < 2 && !is_boundary; ++k) {
                                        coord_type x = abc;

                                        //x[bidx0] -= j;
                                        //x[bidx1] -= k;
                                        x[d] += 2*k-1;
                                        if(!g(x) && og(x)) {
                                            is_boundary = true;
                                        }
                                    }
                                }
                                if(is_boundary) {
                                    bedges[bidx0][abc[bidx0]].insert(e);
                                    bedges[bidx1][abc[bidx1]].insert(e);
                                }
                            }
                        }
                    }
                }
            }
        }
        return {edges,bedges};
    }
    auto AdaptiveGridFactory::get_edge(const coord_type& start, int jump, int dim) const -> Edge {
        coord_type end = start;
        end[dim] += jump;
        return Edge{{vertex_index(start),vertex_index(end)}};
    }
    auto AdaptiveGridFactory::make_edge_shapes(const coord_type& coord) const -> std::array<coord_type,3> {

        std::array<coord_type,3> edge_shapes{{coord,coord,coord}};
        edge_shapes[0][1]++;
        edge_shapes[0][2]++;
        edge_shapes[1][0]++;
        edge_shapes[1][2]++;
        edge_shapes[2][0]++;
        edge_shapes[2][1]++;
        return edge_shapes;
    }

    auto AdaptiveGridFactory::empty_edge_masks(int level, bool value) const -> std::array<GridData3,3> {
        coord_type shape = levels[level].shape();
        mtao::eigen::stl2eigen(shape) *= width;
        return empty_edge_masks(shape, value);
    }
    auto AdaptiveGridFactory::empty_edge_masks(const coord_type& shape, bool value) const -> std::array<GridData3,3> {
        auto EMS = make_edge_shapes(shape);
        std::array<GridData3,3> edge_masks {{
            GridData3::Constant(value, EMS[0]),
                GridData3::Constant(value, EMS[1]),
                GridData3::Constant(value, EMS[2])}};
        return edge_masks;
    }

    auto AdaptiveGridFactory::make_edge_mask(const GridData3b& mask, int level) const -> std::array<GridData3,3>{

        auto edge_masks = empty_edge_masks(level);
        for(int a = 0; a < mask.shape()[0]; ++a) {
            for(int b = 0; b < mask.shape()[1]; ++b) {
                for(int c = 0; c < mask.shape()[2]; ++c) {
                    auto&& m = mask(a,b,c);
                    if(m.any()) {
                        coord_type abc{{width*a,width*b,width*c}};
                        for(int k = 0; k < 3; ++k) {
                            for(int a = abc[0]; a < (abc[0] + width + ((k == 0)?0:1)); ++a) {
                                for(int b = abc[1]; b < (abc[1] + width + ((k == 1)?0:1)); ++b) {
                                    for(int c = abc[2]; c < (abc[2] + width + ((k == 2)?0:1)); ++c) {

                                        edge_masks[k](a,b,c) = true;
                                    }
                                }
                            }
                        }
                        set_edge_masks(cmask_indexer.shape(),abc, edge_masks, [&](const coord_type& abc) {
                                return m[cmask_index(abc)];
                                });
                    }
                }
            }
        }
        return edge_masks;
    }
    auto AdaptiveGridFactory::make_edge_mask(const GridData3& mask, int level) const -> std::array<GridData3,3>{

        auto edge_masks = empty_edge_masks(mask.shape(),true);

        set_edge_masks(mask.shape(),coord_type{{0,0,0}}, edge_masks,[&mask](const coord_type& abc) {
                return mask(abc);
                });
        std::cout << level << "================" << std::endl;
        print_grid(mask);
        for(auto&& em: edge_masks) {
            print_gridb(em);
        }
        std::cout << "================" << std::endl;
        return edge_masks;
    }

    auto AdaptiveGridFactory::make_edges(const GridData3b& mask, int level) const -> std::tuple<std::array<std::set<Edge>,3>,AxialBEdgeMap> {
        return get_edges(make_edge_mask(mask,level),level);
    }
    auto AdaptiveGridFactory::make_edges(const GridData3& mask, int level) const ->  std::tuple<std::array<std::set<Edge>,3>,AxialBEdgeMap> {
        return get_edges(make_edge_mask(mask,level),level);
    }




    AdaptiveGridFactory::AdaptiveGridFactory(const GridData3& mask): original(!mask) {
        std::array<int,3> shape = mask.shape();
        int size = *std::max_element(shape.begin(),shape.end());
        int level_count = int(std::ceil(std::log2(size)));
        size = 1 << level_count ;
        level_count /= logwidth;
        levels.resize(level_count );
        levels_mask.resize(level_count );
        levels[0] = GridData3b::Constant(0,size/width,size/width,size/width);
        levels_mask[0] = original;
        int init_stride = size / width;
        for(int i = 0; i < init_stride; ++i) {
            for(int j = 0; j < init_stride; ++j) {
                for(int k = 0; k < init_stride; ++k) {
                    for(int ii = 0; ii < width; ++ii) {
                        int I = width * i + ii;
                        if(I >= mask.shape()[0]) {
                            break;
                        }
                        for(int jj = 0; jj < width; ++jj) {
                            int J = width * j+jj;
                            if(J >= mask.shape()[1]) {
                                break;
                            }
                            for(int kk = 0; kk < width; ++kk) {
                                int K = width * k+kk;
                                if(K >= mask.shape()[2]) {
                                    break;
                                }
                                levels[0](i,j,k)[cmask_index(ii,jj,kk)] = original(I,J,K);

                            }
                        }
                    }
                }
            }
        }
        for(auto&& s: shape) {
            s = size+1;
        }

        vertex_indexer = Indexer(shape);
        for(int i = 0; i < levels.size()-1; ++i) {
            auto& me = levels[i];
            auto& next = levels[i+1];
            auto& next_mask = levels_mask[i+1];
            std::array<int,3> shape = me.shape();
            next_mask = GridData3::Constant(false,shape);
            for(auto&& s: shape) {
                s /= width;
            }
            next = GridData3b::Constant(0,shape);



            for(int aa = 0; aa < shape[0]; ++aa) {
                for(int bb = 0; bb < shape[1]; ++bb) {
                    for(int cc = 0; cc < shape[2]; ++cc) {
                        ActiveMask v = 0;
                        for(int a = 0; a < width; ++a) {
                            int aaa = width * aa + a;
                            for(int b = 0; b < width; ++b) {
                                int bbb = width * bb + b;
                                for(int c = 0; c < width; ++c) {
                                    int ccc = width * cc + c;
                                    next_mask(aaa,bbb,ccc) = v[cmask_indexer.index(a,b,c)] = me(aaa,bbb,ccc).any();
                                }
                            }
                        }
                        next(aa,bb,cc) = v;
                    }
                }
            }
        }


    }
    void AdaptiveGridFactory::make_cells(const std::optional<int>& max_level) {

        cells.clear();

        for(auto&& [level,g]: mtao::iterator::enumerate(levels)) {

            if(max_level && *max_level == level) {
                if(level == 0) {

                    make_cells( original, 0);
                    break;
                } else {
                    auto& g = levels[level-1];
                    GridData3 A(g.shape());
                    std::transform(g.begin(),g.end(),A.begin(),[&](const ActiveMask& mask) -> bool {
                            return !mask.none();
                            });
                    make_cells( A, level);
                    break;
                }
            } else {

                for(int a = 0; a < g.shape()[0]; ++a) {
                    for(int b = 0; b < g.shape()[1]; ++b) {
                        for(int c = 0; c < g.shape()[2]; ++c) {
                            auto&& m = g(a,b,c);
                            if(m.any()) {
                                make_cells( g(a,b,c), level, coord_type{{a,b,c}});
                            }
                        }
                    }
                }
            }

        }
        auto cell_grid = grid_from_cells(cells);
        print_grid(cell_grid);
    }
    auto AdaptiveGridFactory::grid_from_cells(const std::map<int,Cell>& cells) const -> GridData3i {

        return AdaptiveGrid::grid_from_cells(original.shape(),cells);
    }

    void AdaptiveGridFactory::add_cell(const coord_type& abc, int jump) {
        int cid = original.index(abc);
        assert(cells.find(cid) == cells.end());
        cells[cid] = std::make_tuple(abc,jump);
    }

    void AdaptiveGridFactory::make_cells(const ActiveMask& mask, int level, const coord_type& coord) {
        int jump = 1 << (level*logwidth);
        int pjump = 1 << ((level+1)*logwidth);
        mtao::Vec3i base = pjump * mtao::eigen::stl2eigen(coord);

        for(int a = 0; a < width; ++a) {
            for(int b = 0; b < width; ++b) {
                for(int c = 0; c < width; ++c) {
                    coord_type abc{{a,b,c}};
                    if(!mask[cmask_index(abc)]) {
                        mtao::eigen::stl2eigen(abc) = jump * mtao::eigen::stl2eigen(abc) + base;
                        add_cell(abc,jump);
                    }
                }
            }
        }
    }
    void AdaptiveGridFactory::make_cells(const GridData3& mask, int level) {
        int jump = 1 << (level*logwidth);

        for(int a = 0; a < mask.shape()[0]; ++a) {
            for(int b = 0; b < mask.shape()[1]; ++b) {
                for(int c = 0; c < mask.shape()[2]; ++c) {
                    coord_type abc{{a,b,c}};
                    if(!mask(abc)) {
                        mtao::eigen::stl2eigen(abc) = jump * mtao::eigen::stl2eigen(abc);
                        add_cell(abc,jump);
                    }
                }
            }
        }
    }
    auto AdaptiveGrid::grid() const -> GridData3i {
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
                            std::cout << "Overlapping cells!" << std::endl;
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

    AdaptiveGrid AdaptiveGridFactory::create() const {
        return AdaptiveGrid(original.shape(),cells);

    }
    auto AdaptiveGrid::boundary_triplets(int offset) const -> std::vector<Eigen::Triplet<double>> {
        std::vector<Eigen::Triplet<double>> trips;
        auto& b = boundary();
        for(auto&& [i,e]: mtao::iterator::enumerate(b)) {
            if(is_valid_edge(e)) {
                auto [l,h] = e;
                trips.emplace_back(offset+i,l,-1);
                trips.emplace_back(offset+i,h,1);
            }
        }
        return trips;
    }
    auto AdaptiveGrid::boundary(const GridData3i& grid) const -> std::vector<Edge> {
        std::set<Edge> boundaries;
        auto add_bdry = [&](int a, int b) {
            Edge e{{a,b}};
            if(is_valid_edge(e)) {
                boundaries.emplace(e);
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
                                add_bdry(-1,grid(abc));
                            } else {
                                coord_type tmp = abc;
                                tmp[d]--;
                                add_bdry(grid(tmp),-1);
                            }
                        } else {
                            int pidx = cell_index(abc);
                            coord_type tmp = abc;
                            tmp[d]--;
                            int nidx = cell_index(tmp);
                            int pid = grid.get(pidx);
                            int nid = grid.get(nidx);
                            if(pid >= 0 && nid >= 0) {
                                add_bdry(nid,pid);
                            }
                        }
                    }
                }
            }
        }

        std::vector<Edge> boundaries_vec(boundaries.size());
        std::copy(boundaries.begin(),boundaries.end(),boundaries_vec.begin());

        return boundaries_vec;
    }
    void AdaptiveGrid::make_boundary() {
        m_boundary = boundary(grid());
    }

    mtao::VecXd AdaptiveGrid::dual_edge_lengths() const {

        auto&dx = Base::dx();
        mtao::VecXd ret(m_boundary.size());
        for(auto&& [i,e]: mtao::iterator::enumerate(m_boundary)) {
            if(!is_valid_edge(e)) continue;
            auto [a,b] = e;
            ret(i) = (dx.asDiagonal() * (cell(a).center() - cell(b).center())).norm();
        }
        return ret;
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
        mtao::VecXd ret(m_boundary.size());
        for(auto&& [i,e]: mtao::iterator::enumerate(m_boundary)) {
            if(!is_valid_edge(e)) continue;
            auto [a,b] = e;
            int k = 0;

            auto& ca = cell(a);
            auto& cb = cell(b);
            mtao::Vec3d cd = (ca.center() - cb.center()).cwiseAbs();
            int minwidth = std::min(ca.width(),cb.width());
            cd.maxCoeff(&k);
            double v = minwidth * minwidth;
            for(int j = 0; j < 2; ++j) {
                v *= dx((j+k)%3);
            }
            ret(i) = v;
        }
        return ret;
    }
    int AdaptiveGrid::num_faces() const {
        return boundary().size();
    }
    int AdaptiveGrid::num_cells() const {
        return cells().size();
    }
}
