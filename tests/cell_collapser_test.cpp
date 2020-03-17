#include <iostream>
#include <mandoline/construction/cell_collapser.hpp>
#include <mtao/linear_algebra/gram_schmidt.hpp>
#include <mtao/geometry/mesh/face_normals.hpp>
#include <catch2/catch.hpp>
#include <spdlog/spdlog.h>
#include <iterator>
#include "debug_cutface.hpp"

using E = std::array<int, 2>;
using namespace mandoline::construction;
using namespace mandoline;


TEST_CASE("TetOrientation", "[cell_collapser]") {
    mtao::ColVecs3d V(3,4);
    V.col(0) << 0,0,0;
    V.col(1) << 1,0,0;
    V.col(2) << 0,1,0;
    V.col(3) << 0,0,1;
    mtao::ColVecs3i F(3,4);
    F.col(0) << 0,2,3;
    F.col(1) << 0,3,1;
    F.col(2) << 0,1,2;
    F.col(3) << 1,2,3;

    mtao::ColVecs3d N(3,4);
    N.col(0) << 1,0,0;
    N.col(1) << 0,1,0;
    N.col(2) << 0,0,1;
    N.col(3) << 1,1,1;

    N.col(3).normalize();
    auto NN = mtao::geometry::mesh::face_normals(V,F);
    std::cout << N << std::endl << std::endl;
    std::cout << NN << std::endl;
    std::map<int,CutFace<3>> faces;
    for(int i = 0; i < 4; ++i) {

        auto f = F.col(i);
        faces[i] = CutFace<3>(
                coord_mask<3>{},
                {f(0),f(1),f(2)},// 0 1 2 3
                i,
                N.col(i)
                );
    }
    for(auto&& [idx,f]: faces) {
        CHECK(check_convex_face_normal(V,f));
    }

    CellCollapser cc(faces);


    cc.merge(V,false);

}

TEST_CASE("Windmill", "[cell_collapser]") {


    for(int N = 2; N < 4; ++N) {
        mtao::ColVecs3d V(3,N+2);
        std::map<int,CutFace<3>> faces;
        int face_count = 0;
        V.col(0) = mtao::Vec3d(0,0,-1);
        V.col(1) = mtao::Vec3d(0,0,1);
        mtao::VecXd ang = mtao::VecXd::LinSpaced(N, 0, 2 * M_PI * (N - 1.) / N);
        V.row(0).tail(N) = ang.array().cos().transpose();
        V.row(1).tail(N) = ang.array().sin().transpose();
        V.row(2).tail(N).setZero();
        std::cout << V << std::endl;


        // random frame ( chose GS to generate it because i wanted to write one, not becuase it's a ogod choice )
        //mtao::Mat3d R;
        //R.setRandom();
        //mtao::linear_algebra::gram_schmidt_in_place(R);
        //V = R * V;

        auto make_face = [&](int ) {
        };
        for(int i = 0; i < N; ++i) {
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>{},
                    {0,1,i+2},// 0 1 2 3
                    face_count,
                    (V.col(1)-V.col(0)).cross(V.col(i+2)-V.col(0)).normalized()
                    );
        }
        for(auto&& [idx,f]: faces) {
            CHECK(check_convex_face_normal(V,f));
        }

        CellCollapser cc(faces);


        cc.merge(V,false);
        //for(auto&& [hf,bd]: cc.m_halfface_to_cell) {
        //    auto&& [face_index, edge] = hf;
        //    auto&& [a,b] = edge;
        //    auto&& [bs,sgn] = bd;
        //    spdlog::warn("{}: ({},{}) => {}({})", face_index,a,b,bs,sgn);
        //}


        //for(auto&& [idx,mp]: mtao::iterator::enumerate(cc.cell_boundaries())) {
        //    std::cout << "Cell " << idx << ")";
        //    for(auto&& [fidx,sgn]: mp) {
        //        std::cout << (sgn?'-':'+') << fidx << " ";
        //    }
        //    std::cout << std::endl;
        //}
        auto cf = cc.cell_faces();
        for(auto&& c: cf) {
            REQUIRE(c.size() == 2);
            auto it = c.begin();
            int a = *it++;
            int b = *it;
            if(b == N-1) {
                if(a != 0){
                    REQUIRE(a == N-2 );
                }
            } else {
                REQUIRE( a+1==b );
            }
        }
    }

}

TEST_CASE("WindmillClosed", "[cell_collapser]") {


    for(int N = 3; N < 10; ++N) {
        mtao::ColVecs3d V(3,N+2);
        std::map<int,CutFace<3>> faces;
        int face_count = 0;
        V.col(0) = mtao::Vec3d(0,0,-1);
        V.col(1) = mtao::Vec3d(0,0,1);
        mtao::VecXd ang = mtao::VecXd::LinSpaced(N, 0, 2 * M_PI * (N - 1.) / N);
        V.row(0).tail(N) = ang.array().cos().transpose();
        V.row(1).tail(N) = ang.array().sin().transpose();
        V.row(2).tail(N).setZero();


        // random frame ( chose GS to generate it because i wanted to write one, not becuase it's a ogod choice )
        mtao::Mat3d R;
        R.setRandom();
        mtao::linear_algebra::gram_schmidt_in_place(R);
        //V = R * V;

        auto make_face = [&](int a, int b) {
            int x = a, y = b+2, z = (b+1)%N+2;
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>{},
                    {x,y,z},
                    face_count,
                    (V.col(y)-V.col(x)).cross(V.col(z)-V.col(x))
                    );
        };
        for(int i = 0; i < N; ++i) {
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>{},
                    {0,1,i+2},// 0 1 2 3
                    face_count,
                    (V.col(1)-V.col(0)).cross(V.col(i+2)-V.col(0))
                    );
            make_face(0,i);
            make_face(1,i);
        }

        //for(int i = 0; i < V.cols(); ++i) {
        //    std::cerr << "v " << V.col(i).transpose() << std::endl;
        //}
        //for(auto&& [idx,f]: faces) {
        //    std::cerr << "f ";
        //    for(auto&& c: f.indices) for(auto&& i: c) {
        //        std::cerr << i+1 << " ";
        //    }
        //    std::cerr << std::endl;
        //}

        for(auto&& [idx,f]: faces) {
            CHECK(check_convex_face_normal(V,f));
        }
        CellCollapser cc(faces);


        cc.merge(V,false);
        //for(auto&& [hf,bd]: cc.m_halfface_to_cell) {
        //    auto&& [face_index, edge] = hf;
        //    auto&& [a,b] = edge;
        //    auto&& [bs,sgn] = bd;
        //    fmt::print("{}: ({},{}) => {}({})\n", face_index,a,b,bs,sgn);
        //}


        auto cb = cc.cell_boundaries();
        REQUIRE(cb.size() == N+1);

        for(auto&& [idx,mp]: mtao::iterator::enumerate(cb)) {
            if(mp.size() == 2 * N) {
                // outside negative area cell: the external edges are paired next to one anotehr
                auto it = mp.begin();
                for(int k = 0; k < N; ++k) {
                    int a = std::get<0>(*it++);
                    int b = std::get<0>(*it++);
                    REQUIRE(a == 3*k+1);
                    REQUIRE(b == 3*k+2);
                }
            } else {
                // interior cells are 4 consecutive cells
                int end = std::get<0>(*mp.rbegin());
                const bool last_one = end == 3*(N)-1;
                if(!last_one) {
                    REQUIRE((end % 3) == 0);
                }
                for(auto&& [i,v]: mtao::iterator::enumerate(mp)) {
                    if(last_one && i == 0) {
                        REQUIRE(std::get<0>(v) == 0);
                    } else {
                        REQUIRE(end - 3 + i == std::get<0>(v));
                        //}
                }
            }
        }
    }
}

}
TEST_CASE("TwoCones", "[cell_collapser]") {


    for(int N = 3; N < 10; ++N) {
        mtao::ColVecs3d V(3,N+2);
        std::map<int,CutFace<3>> faces;
        int face_count = 0;
        V.col(0) = mtao::Vec3d(0,0,-1);
        V.col(1) = mtao::Vec3d(0,0,1);
        mtao::VecXd ang = mtao::VecXd::LinSpaced(N, 0, 2 * M_PI * (N - 1.) / N);
        V.row(0).tail(N) = ang.array().cos().transpose();
        V.row(1).tail(N) = ang.array().sin().transpose();
        V.row(2).tail(N).setZero();
        //std::cout << V << std::endl;

        auto make_face = [&](int a, int b, int c) {
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>{},
                    {a,b,c},// 0 1 2 3
                    face_count,
                    (V.col(b)-V.col(a)).cross(V.col(c)-V.col(a)).normalized()
                    );
        };
        for(int i = 0; i < N; ++i) {
            make_face(0,i+2,(i+1)%N+2);
            make_face(1,i+2,(i+1)%N+2);
        }
        std::vector<int> center(N);
        std::iota(center.begin(),center.end(),2);
        faces[face_count++] = CutFace<3>(
                coord_mask<3>{},
                center,
                face_count,
                mtao::Vec3d::UnitZ()
                );

        for(auto&& [idx,f]: faces) {
            CHECK(check_convex_face_normal(V,f));
        }
        CellCollapser cc(faces);


        //std::cout << std::endl;
        cc.bake(V);
        //for(auto&& [hf,bd]: cc.m_halfface_to_cell) {
        //    auto&& [face_index, edge] = hf;
        //    auto&& [a,b] = edge;
        //    auto&& [bs,sgn] = bd;
        //    spdlog::warn("{}: ({},{}) => {}({})", face_index,a,b,bs,sgn);
        //}


        auto cf = cc.cell_faces();
        REQUIRE(cf.size() == 3);
        //std::cout << "num cells: " << cf.size() << std::endl;
        for(auto&& c: cf) {
            if(c.size() == 2*N) {
                for(auto&& [idx,v]: mtao::iterator::enumerate(c)) {
                    REQUIRE(idx==v);
                }
            } else {
                int sgn = *c.begin() % 2;
                REQUIRE(c.size() == N+1);
                for(auto&& c: c) {
                    if(c != 2 * N) {
                        REQUIRE(sgn == c%2);
                    }
                }
            }
            //std::copy(c.begin(),c.end(),std::ostream_iterator<int>(std::cout,","));
            //std::cout << std::endl;
        }
    }

}

TEST_CASE("Cubes", "[cell_collapser]") {


    for (int N = 2; N < 5; ++N) {
        std::stringstream ss;
        ss << "Cubes: " << N;
        SECTION(ss.str());
        spdlog::info("{} Cubes in a line", N);

        mtao::ColVecs3d V(3,4*N+4);
        std::map<int,CutFace<3>> faces;
        int face_count = 0;
        for(int i = 0; i <= N; ++i) {
            V.col(4*i+0) = mtao::Vec3d(0,0,i);
            V.col(4*i+1) = mtao::Vec3d(1,0,i);
            V.col(4*i+2) = mtao::Vec3d(1,1,i);
            V.col(4*i+3) = mtao::Vec3d(0,1,i);
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>(2,i),
                    {4*i,4*i+1,4*i+2,4*i+3},// 0 1 2 3
                    std::array<int,2>{{2,i}},
                    mtao::Vec3d::UnitX().cross(mtao::Vec3d::UnitY())

                    );
        }
        for(int i = 0; i < N; ++i) {
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>(1,0),
                    {4*i+1,4*i+0,4*(i+1)+0,4*(i+1)+1},// 0 1 5 4
                    //{4*i+1,4*i+0,4*(i+1)+0,4*(i+1)+1},// 1 0 4 5
                    std::array<int,2>{{1,0}},
                    mtao::Vec3d::UnitZ().cross(mtao::Vec3d::UnitX())
                    );
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>(1,1),
                    {4*i+2,4*i+3,4*(i+1)+3,4*(i+1)+2},// 3 2 6 7
                    //{4*i+2,4*i+3,4*(i+1)+3,4*(i+1)+2},// 2 3 7 6
                    std::array<int,2>{{1,1}},
                    mtao::Vec3d::UnitZ().cross(mtao::Vec3d::UnitX())
                    );
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>(0,0),
                    {4*i+0,4*i+3,4*(i+1)+3,4*(i+1)+0},// 0 3 7 4
                    std::array<int,2>{{0,0}},
                    mtao::Vec3d::UnitY().cross(mtao::Vec3d::UnitZ())
                    );
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>(0,1),
                    {4*i+1,4*i+2,4*(i+1)+2,4*(i+1)+1},// 1 2 6 5
                    std::array<int,2>{{0,1}},
                    mtao::Vec3d::UnitY().cross(mtao::Vec3d::UnitZ())
                    );
        }

        for(auto&& [idx,f]: faces) {
            auto&& indices = *f.indices.begin();
            mtao::ColVecs2d v(2,indices.size());
            int d = f.as_axial_axis();
            for(int j = 0; j < indices.size(); ++j) {
                v(0,j) = V((d+1)%3,indices[j]);
                v(1,j) = V((d+2)%3,indices[j]);
            }
            /*
               for(int j = 0; j < indices.size(); ++j) {
               auto a = V.col(indices[(j+0)%indices.size()]);
               auto b = V.col(indices[(j+1)%indices.size()]);
               auto c = V.col(indices[(j+2)%indices.size()]);
               mtao::Vec3d n = (b-a).cross(c-a);
               for(int k = 0; k < 3; ++k) {
               REQUIRE(f.N(k) == Approx(n(k)));
               }
               }
               */
            double vol = mtao::geometry::curve_volume(v);
            REQUIRE( vol == Approx(1.)); // double checking the test has the right orientations
            //REQUIRE( std::signbit(vol) == ((d%2==1) ^ std::signbit(f.N(d)))); // double checking the test has the right orientations
            //std::cout << idx << ")" << d << ": " << vol << std::endl;
        }

        //std::cout << std::endl;

        for(auto&& [idx,f]: faces) {
            CHECK(check_convex_face_normal(V,f));
        }
        CellCollapser cc(faces);


        //std::cout << std::endl;
        cc.bake(V);
        //for(auto&& [hf,bd]: cc.m_halfface_to_cell) {
        //    auto&& [face_index, edge] = hf;
        //    auto&& [a,b] = edge;
        //    auto&& [bs,sgn] = bd;
        //    spdlog::warn("{}: ({},{}) => {}({})", face_index,a,b,bs,sgn);
        //}


        auto cf = cc.cell_faces();
        std::cout << "num cells: " << cf.size() << std::endl;
        REQUIRE(cf.size() == N+1);
        for(auto&& c: cf) {
            /*
            if(c.size() == 2*N) {
                for(auto&& [idx,v]: mtao::iterator::enumerate(c)) {
                    REQUIRE(idx==v);
                }
            } else {
                int sgn = *c.begin() % 2;
                REQUIRE(c.size() == N+1);
                for(auto&& c: c) {
                    if(c != 2 * N) {
                        REQUIRE(sgn == c%2);
                    }
                }
            }
            */
            std::copy(c.begin(),c.end(),std::ostream_iterator<int>(std::cout,","));
            std::cout << std::endl;
        }

    }

}

