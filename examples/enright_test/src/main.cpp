#include "mtao/opengl/Window.h"
#include <iostream>
#include "imgui.h"
#include <memory>
#include <algorithm>
#include "mtao/geometry/mesh/boundary_facets.h"
#include "mtao/geometry/mesh/sphere.hpp"
#include "mtao/geometry/mesh/read_obj.hpp"
#include "mtao/geometry/bounding_box.hpp"
#include "mtao/opengl/drawables.h"
#include <mtao/types.h>
#include "enright.h"
#include "mtao/geometry/mesh/eltopo.h"
#include "mtao/geometry/mesh/lostopos.hpp"

#include <Corrade/Utility/Arguments.h>
#include "mtao/opengl/objects/mesh.h"
#include <mandoline/construction/construct_imgui.hpp>
#include <mandoline/construction/remesh_self_intersections.hpp>
#include <mandoline/mesh3.hpp>

#include <glm/gtc/matrix_transform.hpp> 
using TrackerType = ElTopoTracker;
//using TrackerType = mtao::geometry::mesh::LosToposTracker;



class MeshViewer: public mtao::opengl::Window3 {
    public:
        using ColVectors3d = mtao::ColVectors<double,3>;
        using ColVectors3i = mtao::ColVectors<int,3>;
        ColVectors3d V;
        ColVectors3i F;
        std::unique_ptr<TrackerType> tracker;

        //Cutcell mesh parameters
        Eigen::AlignedBox<float,3> bbox;
        std::array<int,3> N{{5,5,5}};
        int& NI=N[0];
        int& NJ=N[1];
        int& NK=N[2];
        bool use_cube = false;
        int adaptive_level=0;

        //Construction helper
        std::optional<mandoline::construction::DeformingGeometryConstructor> constructor;
        //the mesh we're using
        std::optional<mandoline::CutCellMesh<3>> ccm;



        MeshViewer(const Arguments& args): Window3(args), _wireframe_shader{Magnum::Shaders::MeshVisualizer::Flag::Wireframe} {
            Corrade::Utility::Arguments myargs;
            std::tie(V,F) = mtao::geometry::mesh::sphere<double>(4);
            mesh.setTriangleBuffer(V.cast<float>(),F.cast<unsigned int>());
            auto E = mtao::geometry::mesh::boundary_facets(F);

            tracker = std::make_unique<TrackerType>(V,F);



            phong_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{mesh,_shader, drawables()};
            mv_drawable = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{mesh,_wireframe_shader, drawables()};

            mesh.setParent(&root());

        }
        void gui() override {
            if(mv_drawable) {
                mv_drawable->gui();
            }
            if(phong_drawable) {
                phong_drawable->gui();
            }

            if(mv_drawable) {
                mv_drawable->gui();
            }

        }


        void animate() {
            static double t = 0;
            static double dt = 0.01;

            double curt = 0;
            double subdt = dt;
            while(curt < dt) {
                double fine_t = t + curt;
                subdt = dt - curt;

                tracker->improve();


                V = tracker->get_vertices();

                ColVectors3d Vnew = V + (subdt) * 2 * enright_velocities(V, fine_t);

                subdt = tracker->integrate(Vnew,subdt);
                V = Vnew;

                curt = curt + subdt;
                if(curt > dt) {
                    curt = dt;
                }
            }
            auto [V,F] = tracker->get_mesh();
            mesh.setTriangleBuffer(V.cast<float>(),F.cast<unsigned int>());
            t += dt;
            while(t > 1)  {
                t -= 1;
            }
        }

        void draw() override {
            Window3::draw();
            animate();
        }
    private:
        Magnum::Shaders::Phong _shader;
        Magnum::Shaders::MeshVisualizer _wireframe_shader;
        Magnum::Shaders::Flat3D _flat_shader;
        mtao::opengl::objects::Mesh<3> mesh;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* phong_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>* mv_drawable = nullptr;


};




MAGNUM_APPLICATION_MAIN(MeshViewer)
