#include <mtao/types.hpp>
#include <mtao/cmdline_parser.hpp>
#include "mtao/opengl/Window.h"
#include <iostream>
#include "imgui.h"
#include <mtao/eigen/stl2eigen.hpp>
#include "mtao/opengl/shader.h"
#include "mtao/opengl/renderers/bbox.h"
#include "mtao/opengl/renderers/axis.h"
#include <mtao/eigen/stack.h>
#include <memory>
#include <algorithm>

#include <glm/gtc/matrix_transform.hpp> 
#include <glm/gtc/type_ptr.hpp> 
#include "mtao/opengl/renderers/mesh.h"
#include "mtao/opengl/camera.hpp"
#include "mandoline/mesh3.hpp"
using namespace mtao::logging;

using namespace mtao::opengl;

glm::mat4 mvp_it;
bool animate = false;
using Mat = mtao::MatrixX<GLfloat>;
int current_ccm_cell = 0;
int current_ccm_face = 0;
std::optional<std::string> filename;

mandoline::CutCellMesh<3> ccm;


std::unique_ptr<renderers::MeshRenderer> mesh_renderer;
std::unique_ptr<renderers::MeshRenderer> cell_renderer;
std::unique_ptr<renderers::MeshRenderer> face_renderer;

std::unique_ptr<Window> window;


Camera3D cam;

void set_mvp(int w, int h) {
    cam.set_shape(w,h);
    cam.pan();
    cam.zNear() = 0.00001;
    cam.set_perspective();

    cam.update();
}



void update_face() {
    if(ccm.faces().empty()) {
        return;
    }
    current_ccm_face = std::clamp<int>(current_ccm_face,0,ccm.faces().size()-1);
    auto [V,F] = ccm.compact_triangulated_face(current_ccm_face);


    face_renderer->setMesh(V.cast<float>(), F.cast<unsigned int>());
}

void update_cell() {
    if(ccm.cells().empty()) {
        return;
    }
    current_ccm_cell = std::clamp<int>(current_ccm_cell,0,ccm.cell_size()-1);

    auto [V,F] = ccm.compact_triangulated_cell(current_ccm_cell);


    cell_renderer->setMesh(V.cast<float>(), F.cast<unsigned int>());
}

void read_cutmesh(const std::string& filename) {

    ccm = mandoline::CutCellMesh<3>::from_proto(filename);
    ccm.triangulate_faces();
    auto bbox = ccm.bbox();


    mesh_renderer->setMesh((ccm.origV()).cast<float>(), ccm.origF().cast<unsigned int>());
    cell_renderer->setBuffers();
    face_renderer->setBuffers();
    auto m = ((bbox.min() + bbox.max())/2).eval();
    cam.target_pos() = glm::vec3(m.x(),m.y(),m.z());
    cam.camera_pos() = glm::vec3(m.x(),m.y(),m.z() + 2);
    cam.update();
    std::cout << "cut edge made" << std::endl;
}


void prepare_mesh() {

    if(static bool one_use=true; one_use) { one_use = false;
        mesh_renderer = std::make_unique<renderers::MeshRenderer>(3);
        face_renderer = std::make_unique<renderers::MeshRenderer>(3);
        cell_renderer = std::make_unique<renderers::MeshRenderer>(3);

        mesh_renderer->hide_all();
        mesh_renderer->set_face_style(renderers::MeshRenderer::FaceStyle::Phong);
        mesh_renderer->face_color() = glm::vec4(.2,.2,.2,1);

        cell_renderer->hide_all();
        cell_renderer->set_face_style(renderers::MeshRenderer::FaceStyle::Phong);
        cell_renderer->face_color() = glm::vec4(.2,.2,.4,1);

        face_renderer->hide_all();
        face_renderer->set_face_style(renderers::MeshRenderer::FaceStyle::Phong);
        face_renderer->face_color() = glm::vec4(.2,.4,.2,1);

    }






}






ImVec4 clear_color = ImColor(114, 144, 154);






void do_animation() {
}

void render(int width, int height) {
    // Rendering
    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);

    set_mvp(width,height);


    mesh_renderer->set_mvp(cam.mvp());
    mesh_renderer->set_mvp(cam.mv(),cam.p());
    cell_renderer->set_mvp(cam.mvp());
    cell_renderer->set_mvp(cam.mv(),cam.p());
    face_renderer->set_mvp(cam.mvp());
    face_renderer->set_mvp(cam.mv(),cam.p());

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    face_renderer->render();
    cell_renderer->render();
    mesh_renderer->render();


}



void gui_func() {
    {

        auto&& io = ImGui::GetIO();

        ImGui::Checkbox("animate", &animate);
        mesh_renderer->imgui_interface("Mesh Renderer");
        face_renderer->imgui_interface("Face Renderer");
        cell_renderer->imgui_interface("Cell Renderer");


        ImGui::ColorEdit3("clear color", (float*)&clear_color);
        //auto gpos = grid_mouse_pos();
        //glm::ivec2 gposi;
        //gposi.x = std::max<int>(0,std::min<int>(NI-1,gpos.x));
        //gposi.y = std::max<int>(0,std::min<int>(NJ-1,gpos.y));
        //ImGui::Text("Application average %.3f ms/frame (%.1f FPS) [%.3f,%.3f]", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate,gpos.x,gpos.y);
        //ImGui::Text("[%d,%d]: %.5f", gposi.x,gposi.y,data(gposi.x,gposi.y));


        {
            bool dirty = ImGui::SliderInt("Face position",&current_ccm_face,0,ccm.faces().size()-1);
            if(dirty) {
                update_face();
            }
        }
        {
            bool dirty = ImGui::SliderInt("Cell position",&current_ccm_cell,0,ccm.cell_size()-1);
            if(dirty) {
                update_cell();
            }
        }


    }
}
void set_keys() {
    auto&& h= window->hotkeys();
    //h.add([&]() {
    //        prepare_cutmesh();
    //        },"Clear the curves", GLFW_KEY_SPACE);
}


int main(int argc, char * argv[]) {

    mtao::CommandLineParser clp;
    clp.parse(argc,argv);

    if(clp.args().size() > 0) {
        filename = clp.arg(0);
    } else {
        std::cout << "Need a filename" << std::endl;
        return 1;
    }
    set_opengl_version_hints(4,5);
    window = std::make_unique<Window>();
    set_keys();
    window->set_gui_func(gui_func);
    window->set_render_func(render);
    window->makeCurrent();

    prepare_mesh();
    read_cutmesh(*filename);
    window->run();

    return 0;
}


