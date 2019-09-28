#ifdef USE_OLD_STUFF
#include "mtao/opengl/Window.h"
#include "enright.h"
#include <iostream>
#include <chrono>
#include "imgui.h"
#include <memory>
#include <algorithm>
#include "mtao/opengl/renderers/mesh.h"
#include "mtao/geometry/mesh/sphere.hpp"
#include "mtao/geometry/mesh/read_obj.hpp"
#include "mtao/geometry/bounding_box.hpp"
#include "mtao/opengl/camera.hpp"
#include <mtao/types.h>
#include "enright.h"
#include "mtao/geometry/mesh/eltopo.h"
#include "mtao/geometry/mesh/lostopos.hpp"

#include <glm/gtc/matrix_transform.hpp> 

using namespace mtao::opengl;
using namespace mtao::logging;
Camera3D cam;

using ClockType = std::chrono::steady_clock;


bool animate = false;
bool save_frame = false;
std::unique_ptr<Window> window;
std::unique_ptr<renderers::MeshRenderer> renderer;
ImVec4 clear_color = ImColor(114, 144, 154);
double scale = 1.0;
Eigen::Vector3d center;

void prepare_mesh(const ColVectors3d& V, const ColVectors3i&F) {
    renderer = std::make_unique<renderers::MeshRenderer>(3);


    auto bb = mtao::geometry::bounding_box(V);

    scale = bb.sizes().maxCoeff();
    center = bb.center();

    ColVectors3d VV = (V.colwise()-center) /  scale;

    renderer->setMesh(VV.cast<GLfloat>(),F.cast<GLuint>(),true);

    renderers::MeshRenderer::MatrixXgf C = renderer->computeNormals(V.cast<GLfloat>(),F.cast<GLuint>()).array();
    renderer->setColor(C);
    renderer->set_edge_style(renderers::MeshRenderer::EdgeStyle::Mesh);
}
void set_mvp(int w, int h) {
    cam.set_shape(w,h);

    //cam.v() = glm::lookAt(glm::vec3(1,0,0), glm::vec3(0,0,0), glm::vec3(0,1,0));
    cam.pan();
    cam.update();




}
void animate_func() {

}

void gui_func() {
    {

        ImGui::ColorEdit3("clear color", (float*)&clear_color);
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        ImGui::Checkbox("Animate",&animate);
        ImGui::Checkbox("Store",&save_frame);

        renderer->imgui_interface();
    }
    if(ImGui::Button("Reset  Camera?")) {
        cam.reset();
    }

}

void render(int width, int height) {
    // Rendering
    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);

    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);

    set_mvp(width,height);


    renderer->set_mvp(cam.mvp());
    renderer->set_mvp(cam.mv(),cam.p());

    renderer->render();
    if(save_frame) {
        window->save_frame("frame.png");
        save_frame = false;
    }


    if(animate) {

        animate_func();
    }
}



int main(int argc, char * argv[]) {

    set_opengl_version_hints();
    window = std::make_unique<Window>();
    window->set_gui_func(gui_func);
    window->set_render_func(render);
    window->makeCurrent();

    if(argc < 2) {
        std::cout << "Need an obj input!" << std::endl;
        std::cout << "Loading a sphere mesh instead" << std::endl;
        std::tie(V,F) = mtao::geometry::mesh::sphere<double>(3);
    } else {
        std::tie(V,F) = mtao::geometry::mesh::read_objD(argv[1]);
    }
    V *= .15;
    V.colwise() += Eigen::Vector3d(.35,.35,.35);
    tracker = std::make_unique<TrackerType>(V,F);
    prepare_mesh(V,F);

    animate = true;
    /*
    window->record([&](int frame) -> bool {
            return frame < 20;
            }, "frame");
            */
    window->run();

    exit(EXIT_SUCCESS);
}


#else

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
#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/ArrayView.h>

#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Shaders/Phong.h>
#include <Corrade/Utility/Arguments.h>
#include "mtao/opengl/objects/mesh.h"

#include <glm/gtc/matrix_transform.hpp> 
//using TrackerType = ElTopoTracker;
using TrackerType = mtao::geometry::mesh::LosToposTracker;



class MeshViewer: public mtao::opengl::Window3 {
    public:
        using ColVectors3d = mtao::ColVectors<double,3>;
        using ColVectors3i = mtao::ColVectors<int,3>;
        ColVectors3d V;
        ColVectors3i F;

        std::unique_ptr<TrackerType> tracker;
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
    }


    void animate() {
        static double t = 0;
        static double dt = 0.001;

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
#endif
