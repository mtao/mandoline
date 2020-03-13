#include <iostream>
#include <mandoline/construction/cell_collapser.hpp>
#include <catch2/catch.hpp>
#include <iterator>

using E = std::array<int, 2>;
using namespace mandoline::construction;
using namespace mandoline;

TEST_CASE("Cubes", "[cell_collapser]") {


    for (int N = 1; N < 4; ++N) {
        std::stringstream ss;
        ss << "Cubes: " << N;
        SECTION(ss.str());

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
                    {4*i,4*i+1,4*i+2,4*i+3},
                    std::array<int,2>{{2,i}},
                    mtao::Vec3d::UnitZ()
                    );
        }
        for(int i = 0; i < N; ++i) {
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>(1,0),
                    {4*i+0,4*i+1,4*(i+1)+1,4*(i+1)+0},
                    std::array<int,2>{{1,0}},
                    mtao::Vec3d::UnitY()
                    );
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>(1,1),
                    {4*i+3,4*i+2,4*(i+1)+2,4*(i+1)+3},
                    std::array<int,2>{{1,1}},
                    mtao::Vec3d::UnitY()
                    );
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>(0,0),
                    {4*i+0,4*i+1,4*(i+1)+1,4*(i+1)+0},
                    std::array<int,2>{{0,0}},
                    mtao::Vec3d::UnitX()
                    );
            faces[face_count++] = CutFace<3>(
                    coord_mask<3>(0,1),
                    {4*i+3,4*i+2,4*(i+1)+2,4*(i+1)+3},
                    std::array<int,2>{{0,1}},
                    mtao::Vec3d::UnitX()
                    );
        }

        for(auto&& [idx,f]: faces) {
            auto&& indices = *f.indices.begin();
            mtao::Vec2d v(2,indices.size());
            int d = f.as_axial_axis();
            for(int j = 0; j < indices.size(); ++j) {
                v(0,j) = V((d+1)%3,indices[j]);
                v(1,j) = V((d+2)%3,indices[j]);
            }
            auto vol = mtao::geometry::curve_volume(V);
        }


        CellCollapser cc(faces);
        cc.bake(V);
        ;
        auto cf = cc.cell_faces();
        std::cout << "num cells: " << cf.size() << std::endl;
        for(auto&& c: cf) {
            std::copy(c.begin(),c.end(),std::ostream_iterator<int>(std::cout,","));
            std::cout << std::endl;
        }

    }

}

