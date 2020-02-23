#include <iostream>
#include <mandoline/construction/face_collapser.hpp>
#include <catch2/catch.hpp>

using E = std::array<int,2>;
using namespace mandoline::construction;

TEST_CASE("Polygon", "[face_collapser]") {
    for(int N = 3; N < 50; ++N) {
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
            auto faces = fc.faces_no_holes();
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
            auto faces = fc.faces_no_holes();
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
    for(int N = 3; N < 50; ++N) {
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
            auto faces = fc.faces_no_holes();
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
            auto faces = fc.faces_no_holes();
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
TEST_CASE("Near Axis Fan", "[face_collapser]") {
    for(int N = 3; N < 50; ++N) {
        std::stringstream ss;
        SECTION(ss.str());
        mtao::ColVecs2d V(2,N+1);

        V.row(0).head(N).setConstant(1e-8);
        V.row(1).head(N) = mtao::VecXd::LinSpaced(N,1,1e10);
        V.col(N).setZero();


        std::set<E> edges;
        for(int j = 0; j < N; ++j) {
            if(j+1 < N) {
                edges.emplace(E{{j,(j+1)}});
            }
            edges.emplace(E{{j,N}});
        }


        FaceCollapser fc(edges);
        fc.set_edge_for_removal(E{{1,0}});
        fc.bake(V);

        {
            auto faces = fc.faces_no_holes();
            REQUIRE(faces.size() == N-1);
            for(auto [a,l]: faces) {
                REQUIRE(l.size() == 3);
                std::sort(l.begin(),l.end());
                auto it = l.begin();
                int first = *it++;
                int second= *it++;
                int third = *it++;

                REQUIRE(first+1 == second);
                REQUIRE(third == N);
            }
        }
    }

}

TEST_CASE("Near Diagonal Fan", "[face_collapser]") {
    for(int N = 3; N < 50; ++N) {
        std::stringstream ss;
        SECTION(ss.str());
        mtao::ColVecs2d V(2,N+1);

        V.row(0).head(N) = mtao::VecXd::LinSpaced(N,1,1e10);
        V.row(1).head(N) = V.row(0).head(N).array() - .5;
        V.col(N).setZero();


        std::set<E> edges;
        for(int j = 0; j < N; ++j) {
            if(j+1 < N) {
                edges.emplace(E{{j,(j+1)}});
            }
            edges.emplace(E{{j,N}});
        }


        FaceCollapser fc(edges);
        fc.set_edge_for_removal(E{{1,0}});
        fc.bake(V);

        {
            auto faces = fc.faces_no_holes();
            REQUIRE(faces.size() == N-1);
            for(auto [a,l]: faces) {
                REQUIRE(l.size() == 3);
                std::sort(l.begin(),l.end());
                auto it = l.begin();
                int first = *it++;
                int second= *it++;
                int third = *it++;

                REQUIRE(first+1 == second);
                REQUIRE(third == N);
            }
        }
    }

}

TEST_CASE("Triangle with hole", "[face_collapser]") {
    // triangle
    int N = 3;

    mtao::ColVecs2d V(2, 2 * N);

    // triangle in a triangle
    mtao::VecXd ang = mtao::VecXd::LinSpaced(N,0,2*M_PI * (N-1.) / N);
    V.row(0).head(N) = 3*ang.array().cos().transpose();
    V.row(1).head(N) = 3*ang.array().sin().transpose();
    V.row(0).tail(N) = ang.array().cos().transpose();
    V.row(1).tail(N) = ang.array().sin().transpose();



    std::set<E> edges;
    for(int i = 0; i < N; ++i) {
        int a = i;
        int b = (i+1)%N;
        int c = a+N;
        int d = b+N;
    
        edges.emplace(E{{a,b}});
        edges.emplace(E{{c,d}});
    }


    FaceCollapser fc(edges);
    fc.set_edge_for_removal(E{{1,0}});
    fc.bake(V);

    {
        auto faces = fc.faces();
        REQUIRE(faces.size() == 2);
        int num_single_loops = 0;
        int num_double_loops = 0;
        for(auto [a,l]: faces) {
            if(l.size() == 2) {// outer one
                num_double_loops++;

                for(auto&& ll: l) {
                    std::sort(ll.begin(),ll.end());
                    int min = ll[0];
                    for(auto&& [i,v]: mtao::iterator::enumerate(ll)) {
                        if(min == 0) {
                            REQUIRE(i==v);
                        } else {
                            REQUIRE(i+N==v);

                        }
                    }
                }
            } else if(l.size() == 1) {//inner one
                num_single_loops++;
                auto& ll = l.first();
                std::sort(ll.begin(),ll.end());
                for(auto&& [i,v]: mtao::iterator::enumerate(ll)) {
                    REQUIRE(i==v);
                }
            }

        }
        REQUIRE( num_single_loops == 1);
        REQUIRE( num_double_loops == 1);
    }

}
