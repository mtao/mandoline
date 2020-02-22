#include <iostream>
#include <mandoline/construction/face_collapser.hpp>
#include <catch2/catch.hpp>

using E = std::array<int,2>;
using namespace mandoline::construction;

TEST_CASE("Polygon", "[face_collapser]") {
    for(int N = 3; N < 5; ++N) {
        std::stringstream ss;
        ss << "Polygon size " << N;
        SECTION(ss.str());
        mtao::ColVecs2d V(2,N);

        mtao::VecXd ang = mtao::VecXd::LinSpaced(N,0,2*M_PI * (N-1.) / N);
        V.row(0) = ang.array().cos().transpose();
        V.row(1) = ang.array().sin().transpose();


        std::set<E> edges;
        for(int j = 0; j < N; ++j) {
            edges.emplace(E{{j,(j+1)%N}});
        }


        FaceCollapser fc(edges);
        fc.bake(V);


        {
            auto faces = fc.faces();
            REQUIRE(faces.size() == 2);
            {
                auto prit = faces.find(0);
                REQUIRE(prit != faces.end());
                REQUIRE(prit->second.size() == N);
            }
            {
                auto prit = faces.find(1);
                REQUIRE(prit != faces.end());
                REQUIRE(prit->second.size() == N);
            }
        }

        fc = FaceCollapser(edges);
        fc.set_edge_for_removal(E{{1,0}});
        fc.bake(V);

        {
            auto faces = fc.faces();
            REQUIRE(faces.size() == 1);
            {
                auto prit = faces.find(0);
                REQUIRE(prit != faces.end());
                REQUIRE(prit->second.size() == N);
            }
            for(auto&& [i,v]: mtao::iterator::enumerate(faces[1])) {
                REQUIRE(i==v);
            }
        }
    }

}
TEST_CASE("Polygon with Center", "[face_collapser]") {
    for(int N = 3; N < 6; ++N) {
        std::stringstream ss;
        ss << "Polygon size " << N;
        SECTION(ss.str());
        mtao::ColVecs2d V(2,N+1);

        mtao::VecXd ang = mtao::VecXd::LinSpaced(N,0,2*M_PI * (N-1.) / N);
        V.row(0).head(N) = ang.array().cos().transpose();
        V.row(1).head(N) = ang.array().sin().transpose();
        V.col(N).setZero();


        std::set<E> edges;
        for(int j = 0; j < N; ++j) {
            edges.emplace(E{{j,(j+1)%N}});
            edges.emplace(E{{j,N}});
        }


        FaceCollapser fc(edges);
        fc.bake(V);


        {
            auto faces = fc.faces();
            REQUIRE(faces.size() == N+1);
            {
                int triangle_count = 0;
                int non_tri_poly_count = 0;
                for(auto&& [a,b]: faces) {
                    int size = b.size();
                    if(size == 3) {
                        triangle_count++;
                    } else {
                        non_tri_poly_count++;
                    }
                }
                if(N > 3) {
                REQUIRE(triangle_count == N);
                REQUIRE(non_tri_poly_count == 1);
                } else {
                    REQUIRE(triangle_count == N+1);
                }
            }
        }

        fc = FaceCollapser(edges);
        fc.set_edge_for_removal(E{{1,0}});
        fc.bake(V);

        {
            auto faces = fc.faces();
            REQUIRE(faces.size() == N);
            for(auto [a,l]: faces) {
                REQUIRE(l.size() == 3);
                std::sort(l.begin(),l.end());
                auto it = l.begin();
                int first = *it++;
                int second= *it++;
                int third = *it++;

                if(first == 0 && second != 1) {
                    REQUIRE(second == N-1);
                } else {
                    REQUIRE(first+1 == second);
                }
                REQUIRE(third == N);


            }
        }
    }

}
