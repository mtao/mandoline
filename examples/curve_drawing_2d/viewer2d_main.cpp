/*
#include <mtao/types.hpp>
#include "mtao/opengl/Window.h"
#include <iostream>
#include "imgui.h"
#include "eigen_utils.hpp"
#include "mtao/opengl/shader.h"
#include "mtao/opengl/renderers/bbox.h"
#include "mtao/opengl/renderers/axis.h"
#include <mtao/eigen/stack.h>
#include <mtao/logging/timer.hpp>
#include <memory>
#include <algorithm>

#include <glm/gtc/matrix_transform.hpp> 
#include <glm/gtc/type_ptr.hpp> 
#include "mtao/opengl/renderers/mesh.h"
#include "mtao/opengl/camera.hpp"
#include <mtao/eigen_utils.h>
#include <mtao/geometry/grid/triangulation.hpp>
#include "plcurve2.hpp"
#include "mandoline/construction/generator.hpp"
#include "mandoline/mesh2.hpp"

using namespace mtao::opengl;
using namespace mandoline;

using GridData = mtao::geometry::grid::GridDataD<double,2>;
using GridDatab = mtao::geometry::grid::GridDataD<bool,2>;
//std::tuple<mtao::ColVectors<float,2>,mtao::ColVectors<unsigned int,3>,mtao::ColVectors<unsigned int,2>> mesh(const std::array<int,2>&);

PLCurve2 curve;

glm::mat4 mvp_it;
bool animate = false;
using Mat = mtao::MatrixX<GLfloat>;

int NI=4;
int& NJ = NI;
int& grid_resolution = NI;
CutCellMesh<2> ccm(std::array<int,2>{{NI,NJ}});

std::unique_ptr<renderers::MeshRenderer> cutedge_renderer;
std::unique_ptr<renderers::MeshRenderer> cutmesh_renderer;
std::unique_ptr<renderers::MeshRenderer> cell_renderer;
std::unique_ptr<renderers::MeshRenderer> edge_renderer;
std::unique_ptr<renderers::MeshRenderer> point_renderer;
std::unique_ptr<renderers::MeshRenderer> grid_renderer;
std::unique_ptr<renderers::MeshRenderer> grid_dual_renderer;
std::unique_ptr<renderers::BBoxRenderer2> bbox_renderer;
std::unique_ptr<renderers::AxisRenderer> axis_renderer;

std::unique_ptr<Window> window;


Camera2D cam;

void set_mvp(int w, int h) {
    cam.set_shape(w,h);
    auto&& m = cam.m();
    m = glm::mat4();

    //cam.v() = glm::lookAt(glm::vec3(1,0,0), glm::vec3(0,0,0), glm::vec3(0,1,0));
    cam.pan();
    cam.update();

    mvp_it = cam.mvp_inv_trans();
}





void set_colors(const GridDatab& data) {
    mtao::RowVectorX<GLfloat> col = Eigen::Map<const mtao::VectorX<bool>>(data.data(),data.size()).cast<float>();

    mtao::MatrixX<GLfloat> Col = mtao::eigen::vstack((col.array()>0).select(col,0),-(col.array()<0).select(col,0),mtao::RowVectorX<GLfloat>::Zero(data.size()));

    Col.array() = Col.array().pow(.5);
    grid_dual_renderer->setColor(Col);

}

void update_edges() {

    auto points = curve.points();
    auto V = points.cast<float>().eval();
    if(V.size() == 0) {
        edge_renderer->setBuffers();
        cutedge_renderer->setBuffers();
        cutmesh_renderer->setBuffers();
    } else {
        auto E = curve.edges();
        edge_renderer->setVertices(V);
        edge_renderer->setEdges(E);
        if(E.size() == 0) {
            return;
        }

        auto  stlp = curve.stl_points();


        construction::CutCellGenerator<2> ccg(stlp,ccm.cell_shape());
        ccg.add_boundary_elements(E.cast<int>());
        ccg.bake();

        ccm = ccg.generate();

        set_colors(ccm.active_grid_cell_mask);

        {
            auto V = ccm.vertices();
            auto E = ccm.cut_edges;
            mtao::ColVectors<float,3> ccg_cols_edge(3,E.cols());
            ccg_cols_edge.setZero();
            for(int i = 0; i < ccg_cols_edge.cols(); ++i) {
                ccg_cols_edge(i%3,i) = 1;
            }
            mtao::ColVectors<float,3> ccg_cols(3,V.cols());
            ccg_cols.setZero();
            for(int i = 0; i < ccm.StaggeredGrid::vertex_size(); ++i) {
                ccg_cols.col(i).setConstant(.7);
            }
            //for(int i: ccm.face_vertices) {
            //    ccg_cols(2,i + ccg.vertex_size()) = 1;
            //}
            //for(int i: ccm.edge_vertices) {
            //    ccg_cols(1,i + ccg.vertex_size()) = 1;
            //}
            //for(int i: ccm.original_vertices) {
            //    ccg_cols(0,i + ccg.vertex_size()) = 1;
            //}
            //for(int i: ccm.original_vertices) {
            //    ccg_cols(2,i) = 1;
            //}
            //for(int i: ccm.new_vertices) {
            //    ccg_cols(0,i) = 1;
            //}
            cutedge_renderer->setVertices(V.cast<float>());
            cutedge_renderer->setEdges(E.cast<unsigned int>());
            cutedge_renderer->setColor(ccg_cols_edge);
            auto F = ccm.faces();
            std::cout << "F size: " << F.rows() << "," << F.cols()  << std::endl;
            cutmesh_renderer->setVertices(V.cast<float>());
            cutmesh_renderer->setFaces(F.cast<unsigned int>());
            cutmesh_renderer->setColor(ccg_cols);
            point_renderer->setVertices(ccm.dual_vertices().rightCols(ccm.cell_size() - ccm.StaggeredGrid::cell_size()).cast<float>());
            //point_renderer->setVertices(ccg.compact_vertices().cast<float>().colwise() - mtao::Vec2f(.5,.5));
            //point_renderer->setVertices(ccm.new_vertices.cast<float>().colwise() - mtao::Vec2f(.5,.5));
            cutedge_renderer->set_vertex_style(renderers::MeshRenderer::VertexStyle::Flat);
            cutedge_renderer->set_face_style();
            cutedge_renderer->set_edge_style(renderers::MeshRenderer::EdgeStyle::Mesh);
        }
    }
}

void prepare_mesh() {


}




ImVec4 clear_color = ImColor(114, 144, 154);



glm::vec2 mouse_pos() {
    return cam.mouse_pos(ImGui::GetIO().MousePos);;
}

glm::vec2 grid_space(const glm::vec2& p) {
    return p * glm::vec2(NI,NJ);
}
glm::vec2 grid_mouse_pos() {
    return grid_space(mouse_pos());
}





void render(int width, int height) {
    return;
    // Rendering
    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);

    set_mvp(width,height);


    grid_dual_renderer->set_mvp(cam.mvp());
    grid_dual_renderer->set_mvp(cam.mv(),cam.p());
    grid_renderer->set_mvp(cam.mvp());
    grid_renderer->set_mvp(cam.mv(),cam.p());
    bbox_renderer->set_mvp(cam.mvp());
    //axis_renderer->set_mvp(cam.mvp() * glm::translate(glm::mat4(),glm::vec3(100,100,0)));
    edge_renderer->set_mvp(cam.mvp());
    cutedge_renderer->set_mvp(cam.mvp());
    cutmesh_renderer->set_mvp(cam.mvp());
    point_renderer->set_mvp(cam.mvp());

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_ALWAYS);
    glEnable(GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);
    grid_renderer->render();
    cutmesh_renderer->render();
    grid_dual_renderer->render();
    
    edge_renderer->render();

    cutedge_renderer->render();
    point_renderer->render();
    glDepthFunc(GL_LESS);
    //bbox_renderer->render();
    //axis_renderer->render();



}


void gui_func() {
    return;
    auto&& io = ImGui::GetIO();

    ImGui::Checkbox("animate", &animate);


    grid_renderer->imgui_interface("Grid Renderer");
    grid_dual_renderer->imgui_interface("Dual Grid Renderer");
    edge_renderer->imgui_interface("Edge Renderer");
    cutedge_renderer->imgui_interface("Cutedge Renderer");
    cutmesh_renderer->imgui_interface("Cutmesh Renderer");
    point_renderer->imgui_interface("Point Renderer");


    ImGui::ColorEdit3("clear color", (float*)&clear_color);
    //auto gpos = grid_mouse_pos();
    auto gpos = mouse_pos();
    auto [c,q] = ccm.coord(Eigen::Map<mtao::Vec2f>(glm::value_ptr(gpos)).cast<double>());
    auto cell = ccm.cell_index(Eigen::Map<mtao::Vec2f>(glm::value_ptr(gpos)).cast<double>());
    auto edge = ccm.nearest_edge_index(Eigen::Map<mtao::Vec2f>(glm::value_ptr(gpos)).cast<double>());
    ImGui::Text("Application average %.3f ms/frame (%.1f FPS) [%.3f,%.3f]", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate,gpos.x,gpos.y);
    //ImGui::Text("[%d,%d]: %.5f", gposi.x,gposi.y,data(gposi.x,gposi.y));
    ImGui::Text("[%d,%d],(%.3f,%.3f): cell %d, edge %d", c[0], c[1],q[0],q[1], cell, edge);


    static int cell2 = 0;
    ImGui::InputInt("Face position2",&cell2);
    if(cell2 >=0 && cell2 < ccm.cell_size() ) {

        auto CIV = ccm.cell(cell2);
        //for(auto&& v: CIV) {
        //    std::cout << v << ",";
        //}
        //std::cout << std::endl;
        mtao::ColVectors<int,2> E(2,CIV.size());
        for(int i = 0; i < CIV.size(); ++i) {
            int a = CIV[i];
            int b = CIV[(i+1)%CIV.size()];
            E.col(i) = mtao::Vec2i(a,b);
        }

        cutmesh_renderer->setEdges(E.cast<unsigned int>());
        if(CIV.size() >= 3) {
            auto F = mtao::geometry::mesh::earclipping(ccm.vertices(),CIV);
            cutmesh_renderer->setFaces(F.cast<unsigned int>());
        }
    }
    ImGui::InputInt("Grid Resolution",&grid_resolution);
    if(ImGui::Button("Restart")) {
        prepare_mesh();
    }






}
void set_keys() {
    auto&& h= window->hotkeys();
    h.add([&]() {
            curve.add_point(glm2eigen(mouse_pos() ).cast<double>());
            update_edges();
            },"Add a point to the curve", GLFW_KEY_P);
    h.add([&]() {
            curve.clear();
            update_edges();
            },"Clear the curves", GLFW_KEY_C);
    h.add([&]() {
            curve.toggle_closed();
            update_edges();
            },"Clear the curves", GLFW_KEY_O);
    h.add([&]() {
            if(curve.size() > 0) {
            curve.pop_back();
            update_edges();
            }
            },"Pop the last curve entry", GLFW_KEY_U);

}


int main(int argc, char * argv[]) {

    set_opengl_version_hints(4,5);
    window = std::make_unique<Window>();
    set_keys();
    window->set_gui_func(gui_func);
    window->set_render_func(render);
    window->makeCurrent();
    cam.set_translation(-glm::vec2(.5,.5));

    prepare_mesh();
    window->run();

    exit(EXIT_SUCCESS);
}
*/



#include "mtao/geometry/mesh/boundary_facets.h"
#include "mtao/geometry/bounding_box.hpp"
#include <mtao/types.hpp>
#include "mtao/opengl/Window.h"
#include <iostream>
#include "imgui.h"
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/eigen/stack.h>
#include <memory>
#include <algorithm>
#include <mtao/opengl/drawables.h>
#include <Magnum/EigenIntegration/Integration.h>
#include <mtao/opengl/objects/mesh.h>
#include <mtao/opengl/objects/bbox.h>
#include <mtao/opengl/objects/grid.h>
#include <Corrade/Utility/Arguments.h>
#include "plcurve2.hpp"
#include "mandoline/construction/generator.hpp"
#include <Magnum/GL/Renderer.h>

#include <thread>
#include "mandoline/mesh3.hpp"
using namespace mtao::logging;
using namespace Magnum::Math::Literals;




class MeshViewer: public mtao::opengl::Window2 {
    public:

        std::optional<mandoline::CutCellMesh<2>> ccm;
        bool show_multi = false;
        int index = 0;
        PLCurve2 curve;
        Eigen::AlignedBox<float,2> bbox;
        std::array<int,2> N{{5,5}};


        mtao::opengl::Vector2 cursor;

        bool show_wireframes = false;

        float scale = 1.1;
        float region_center_scale = 0.0;
        bool do_slice = false;
        mtao::ColVectors<double,4> colors;
        mtao::ColVectors<double,2> cachedV;
        mtao::ColVectors<int,2> cachedE;


        mtao::Vec2f origin = mtao::Vec2f::Zero(), direction = mtao::Vec2f::Unit(2);


        int size() const { return colors.cols(); }



        MeshViewer(const Arguments& args): Window2(args) {
            mtao::logging::make_logger().set_level(mtao::logging::Level::Off);
            mtao::logging::make_logger("profiler").set_level(mtao::logging::Level::Off);
            bbox.min().setConstant(-1);
            bbox.max().setConstant(1);
            //Corrade::Utility::Arguments myargs;
            //myargs.addArgument("filename").parse(args.argc,args.argv);
            //std::string filename = myargs.value("filename");

            update_bbox();


            bbox_drawable = mtao::opengl::make_drawable(bbox_mesh, _flat_shader, drawables());
            bbox_mesh.setParent(&root());
            bbox_drawable->deactivate();
            bbox_drawable->activate_edges();

            bbox_drawable->data().color = 0xffffff_rgbf;


            curve_drawable= new mtao::opengl::Drawable<Magnum::Shaders::Flat2D>{curve_mesh,_flat_shader, curve_drawgroup};
            curve_mesh.setParent(&root());
            curve_drawable->deactivate();
            cutcell_drawable= new mtao::opengl::Drawable<Magnum::Shaders::Flat2D>{cutcell_mesh,_flat_shader, background_drawgroup};
            cutcell_mesh.setParent(&root());
            cutcell_drawable->deactivate();


            grid_mesh.setParent(&root());
            grid_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat2D>(grid_mesh, _flat_shader, background_drawgroup);
            grid_drawable->deactivate();


        }

        void update_bbox() {
            bbox_mesh.set_bbox(bbox);

        }
        void draw() override;
        void gui() override;
        void update_edges();
        void update_curve();

        void mouseMoveEvent(MouseMoveEvent& event) override;
        void mousePressEvent(MouseEvent& event) override;
    private:
        Magnum::SceneGraph::DrawableGroup2D background_drawgroup, curve_drawgroup;


        Magnum::Shaders::Flat2D _flat_shader;
        Magnum::Shaders::VertexColor2D vcolor_shader;
        mtao::opengl::objects::Mesh<2> curve_mesh;
        mtao::opengl::objects::Mesh<2> cutcell_mesh;
        mtao::opengl::objects::Grid<2> grid_mesh;
        mtao::opengl::objects::BoundingBox<2> bbox_mesh;
        mtao::opengl::Drawable<Magnum::Shaders::Flat2D>* bbox_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat2D>* grid_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat2D>* curve_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat2D>* cutcell_drawable = nullptr;

};
void MeshViewer::draw() {
    Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::DepthTest);
    Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);
    Magnum::GL::Renderer::setPointSize(10);

    camera().draw(background_drawgroup);
    Window2::draw();
    camera().draw(curve_drawgroup);
}

void MeshViewer::mouseMoveEvent(MouseMoveEvent& event) {
    Window2::mouseMoveEvent(event);
    cursor = localPosition(event.position());
}
void MeshViewer::mousePressEvent(MouseEvent& event) {
    Window2::mousePressEvent(event);
    if(!ImGui::GetIO().WantCaptureMouse) {
        if(event.button() == MouseEvent::Button::Left) {
            mtao::Vec2d p(cursor.x(),cursor.y());
            curve.add_point(p);
            update_curve();
        }
    }
}

void MeshViewer::gui() {
    if(ImGui::InputInt2("N", N.data()))  {
        //update_bbox();
    }
    if(ImGui::SliderFloat2("min", bbox.min().data(),-2,2))  {
        bbox.min() = (bbox.min().array() < bbox.max().array()).select(bbox.min(),bbox.max());
        update_bbox();
    }
    if(ImGui::SliderFloat2("max", bbox.max().data(),-2,2))  {
        bbox.max() = (bbox.min().array() > bbox.max().array()).select(bbox.min(),bbox.max());
        update_bbox();
    }

    ImGui::Text("Cursor position (%f,%f)", cursor.x(),cursor.y());
}

void MeshViewer::update_edges() {
}
void MeshViewer::update_curve() {

    auto points = curve.points();
    auto V = curve.points().cast<float>().eval();
    auto E = curve.edges();
    if(V.size() == 0) {
        curve_drawable->deactivate();
    } else {
        if(E.size() > 0) {
            curve_mesh.setEdgeBuffer(V,E);
            curve_drawable->activate_edges();
        } else {
            curve_mesh.setVertexBuffer(V);
            curve_drawable->activate_points();
        }
    }

    auto stlp = curve.stl_points();


    auto grid = mtao::geometry::grid::StaggeredGrid2d::from_bbox(bbox.cast<double>(),N);
    mandoline::construction::CutCellGenerator<2> ccg(stlp,grid);
    ccg.add_boundary_elements(E.cast<int>());
    ccg.bake();

    ccm = ccg.generate();

    //set_colors(ccm.active_grid_cell_mask);

    {
        auto V = ccm->vertices();
        auto E = ccm->cut_edges();
        mtao::ColVectors<float,3> ccg_cols_edge(3,E.cols());
        ccg_cols_edge.setZero();
        for(int i = 0; i < ccg_cols_edge.cols(); ++i) {
            ccg_cols_edge(i%3,i) = 1;
        }
        mtao::ColVectors<float,3> ccg_cols(3,V.cols());
        ccg_cols.setZero();
        for(int i = 0; i < ccm->StaggeredGrid::vertex_size(); ++i) {
            ccg_cols.col(i).setConstant(.7);
        }
        //for(int i: ccm.face_vertices) {
        //    ccg_cols(2,i + ccg.vertex_size()) = 1;
        //}
        //for(int i: ccm.edge_vertices) {
        //    ccg_cols(1,i + ccg.vertex_size()) = 1;
        //}
        //for(int i: ccm.original_vertices) {
        //    ccg_cols(0,i + ccg.vertex_size()) = 1;
        //}
        //for(int i: ccm.original_vertices) {
        //    ccg_cols(2,i) = 1;
        //}
        //for(int i: ccm.new_vertices) {
        //    ccg_cols(0,i) = 1;
        //}
        cutcell_mesh.setEdgeBuffer(V.cast<float>().eval(),E.cast<unsigned int>().eval());
        cutcell_drawable->activate_edges();
        //cutedge_renderer->setVertices(V.cast<float>());
        //cutedge_renderer->setEdges(E.cast<unsigned int>());
        //cutedge_renderer->setColor(ccg_cols_edge);
        //auto F = ccm.faces();
        //std::cout << "F size: " << F.rows() << "," << F.cols()  << std::endl;
        //cutmesh_renderer->setVertices(V.cast<float>());
        //cutmesh_renderer->setFaces(F.cast<unsigned int>());
        //cutmesh_renderer->setColor(ccg_cols);
        //point_renderer->setVertices(ccm.dual_vertices().rightCols(ccm.cell_size() - ccm.StaggeredGrid::cell_size()).cast<float>());
        //point_renderer->setVertices(ccg.compact_vertices().cast<float>().colwise() - mtao::Vec2f(.5,.5));
        //point_renderer->setVertices(ccm.new_vertices.cast<float>().colwise() - mtao::Vec2f(.5,.5));
        //cutedge_renderer->set_vertex_style(renderers::MeshRenderer::VertexStyle::Flat);
        //cutedge_renderer->set_face_style();
        //cutedge_renderer->set_edge_style(renderers::MeshRenderer::EdgeStyle::Mesh);

        grid_mesh.set(ccm->vertex_grid());
        grid_drawable->activate_edges();
    }
}



MAGNUM_APPLICATION_MAIN(MeshViewer)

