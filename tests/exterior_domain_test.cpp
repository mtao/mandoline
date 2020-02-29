#include <iostream>
#include <mandoline/exterior_grid.hpp>
#include <catch2/catch.hpp>
#include <iterator>

template <int D>
auto print_eg_dual_edges(const mandoline::ExteriorGrid<D>& eg) {
    auto gb = eg.cell_grid();
    for(auto&& [a,b]: eg.boundary_facet_pairs()) {
        if(a >= 0) {
            auto c = gb.unindex(a);
            std::copy(c.begin(),c.end(),std::ostream_iterator<int>(std::cout," "));
        } else {
            std::cout << a;
        }
        std::cout << " => ";
        if(b >= 0) {
            auto c = gb.unindex(b);
            std::copy(c.begin(),c.end(),std::ostream_iterator<int>(std::cout," "));
        } else {
            std::cout << b;
        }
        std::cout << std::endl;
    }
}

using E = std::array<int, 2>;
using namespace mandoline::construction;

TEST_CASE("2D", "[exterior_grid]") {
    using EG = mandoline::ExteriorGrid<2>;
    using GridDatab = EG::GridDatab;
    
    int N = 5;
    int M = 5;
    GridDatab gb = GridDatab::Constant(true,N,M);

    auto eg = EG(gb);
    REQUIRE(eg.num_faces() == N*(M+1) + (N+1)*M);
    REQUIRE(eg.num_cells() == N*M);

    {
        auto bt = eg.boundary_triplets();
        REQUIRE(bt.size() == 2 * (N-1)*M + 2 * N*(M-1));
        for(auto&& t: bt) {
            REQUIRE(t.row() >= 0);
            REQUIRE(t.row() < eg.num_faces());

            REQUIRE(t.col() >= 0);
            REQUIRE(t.col() < eg.num_cells());
        }

    }
    {
        auto bt = eg.boundary_triplets(true);
        REQUIRE(bt.size() == 4 * eg.num_cells());
        for(auto&& t: bt) {
            REQUIRE(t.row() >= 0);
            REQUIRE(t.row() < eg.num_faces());

            REQUIRE(t.col() >= 0);
            REQUIRE(t.col() < eg.num_cells());
        }

    }
}
TEST_CASE("2D empty", "[exterior_grid]") {
    using EG = mandoline::ExteriorGrid<2>;
    using GridDatab = EG::GridDatab;
    
    int N = 2;
    int M = 2;
    GridDatab gb = GridDatab::Constant(false,N,M);

    auto eg = EG(gb);
    REQUIRE(eg.num_faces() == 0);
    REQUIRE(eg.num_cells() == 0);

    {
        auto bt = eg.boundary_triplets(true);
        REQUIRE(bt.size() == 0);

    }
}
TEST_CASE("2D single ", "[exterior_grid]") {
    using EG = mandoline::ExteriorGrid<2>;
    using GridDatab = EG::GridDatab;
    
    int N = 2;
    int M = 2;
    GridDatab gb = GridDatab::Constant(false,N,M);
    gb(0,0) = 1;

    auto eg = EG(gb);
    print_eg_dual_edges(eg);
    REQUIRE(eg.num_faces() == 2);
    REQUIRE(eg.num_cells() == 1);

    {
        auto bt = eg.boundary_triplets(true);
        REQUIRE(bt.size() == 2);

    }
}
