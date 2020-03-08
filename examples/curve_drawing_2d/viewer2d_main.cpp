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
#include "mandoline/construction/generator2.hpp"
#include <Magnum/GL/Renderer.h>
#include <mtao/geometry/mesh/stack_meshes.hpp>

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

        int face_index = -1;

        mtao::Vec2f origin = mtao::Vec2f::Zero(), direction = mtao::Vec2f::Unit(1);


        int size() const { return colors.cols(); }



        MeshViewer(const Arguments& args): Window2(args) {
            //mtao::logging::make_logger().set_level(mtao::logging::Level::Off);
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
            cutcell_drawable= new mtao::opengl::Drawable<Magnum::Shaders::VertexColor2D>{cutcell_mesh,vcolor_shader, background_drawgroup};
            cutcell_mesh.setParent(&root());
            cutcell_drawable->deactivate();

            cutcell_face_drawable= new mtao::opengl::Drawable<Magnum::Shaders::Flat2D>{cutcell_face_mesh,_flat_shader, background_drawgroup};
            cutcell_face_mesh.setParent(&root());
            cutcell_face_drawable->deactivate();


            grid_mesh.setParent(&root());
            grid_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat2D>(grid_mesh, _flat_shader, background_drawgroup);
            grid_drawable->deactivate();

            reset_curve();

        }

        void update_bbox() {
            bbox_mesh.set_bbox(bbox);

        }
        void draw() override;
        void gui() override;
        void update_edges();
        void update_curve();
        void reset_curve();
        void clear_curve();
        void update_faces();
        void update_ccm();

        void mouseMoveEvent(MouseMoveEvent& event) override;
        void mousePressEvent(MouseEvent& event) override;
        void update_face(int idx);
    private:
        Magnum::SceneGraph::DrawableGroup2D background_drawgroup, curve_drawgroup;


        Magnum::Shaders::Flat2D _flat_shader;
        Magnum::Shaders::VertexColor2D vcolor_shader;
        mtao::opengl::objects::Mesh<2> curve_mesh;
        mtao::opengl::objects::Mesh<2> cutcell_mesh;
        mtao::opengl::objects::Mesh<2> cutcell_face_mesh;
        mtao::opengl::objects::Grid<2> grid_mesh;
        mtao::opengl::objects::BoundingBox<2> bbox_mesh;
        mtao::opengl::Drawable<Magnum::Shaders::Flat2D>* bbox_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat2D>* grid_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat2D>* curve_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::VertexColor2D>* cutcell_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat2D>* cutcell_face_drawable = nullptr;

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

void MeshViewer::clear_curve() {
    curve.clear();
}
void MeshViewer::reset_curve() {

    curve.clear();
    if(!curve.is_closed()) {
        curve.toggle_closed();
    }
    curve.add_point(mtao::Vec2d(-.3,-.3));
    curve.add_point(mtao::Vec2d(-.3,.3));
    curve.add_point(mtao::Vec2d(.3,.3));
    curve.add_point(mtao::Vec2d(.3,-.3));
}

void MeshViewer::gui() {
    if(ImGui::InputInt2("N", N.data()))  {
        //update_bbox();
    }
    if(ImGui::InputInt("Face index", &face_index))  {
        if(ccm) {
            int nf = ccm->num_faces();
            face_index = std::clamp<int>(face_index,0,nf-1);
        }
        update_face(face_index);
    }
    if(ImGui::SliderFloat2("min", bbox.min().data(),-2,2))  {
        bbox.min() = (bbox.min().array() < bbox.max().array()).select(bbox.min(),bbox.max());
        update_bbox();
    }
    if(ImGui::SliderFloat2("max", bbox.max().data(),-2,2))  {
        bbox.max() = (bbox.min().array() > bbox.max().array()).select(bbox.min(),bbox.max());
        update_bbox();
    }

    {
        bool value = curve.is_closed();
        if(ImGui::Checkbox("Closed",&value)) {
            curve.toggle_closed();
            update_curve();
        }
    }
    {
        if(ImGui::Button("Make CCM")) {
            update_ccm();

            update_curve();
        }
        if(ImGui::Button("Clear Curve")) {
            clear_curve();

            update_curve();
        }
        if(ImGui::Button("Reset Curve")) {
            reset_curve();
            update_curve();
        }
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
}

void MeshViewer::update_ccm() {
    auto E = curve.edges();
    auto stlp = curve.stl_points();


    auto grid = mtao::geometry::grid::StaggeredGrid2d::from_bbox(bbox.cast<double>(),N);
    mandoline::construction::CutCellGenerator<2> ccg(stlp,grid);
    ccg.add_boundary_elements(E.cast<int>());
    ccg.bake();

    ccm = ccg.generate();
    std::cout << "HEM data: " << std::endl;
    std::cout << ccm->hem.edges() << std::endl;
    std::cout <<std::endl;

    for(auto&& [a,b]: ccm->exterior_grid.boundary_facet_pairs()) {
        std::cout << a << ":" << b << " ";
    }
    std::cout << std::endl;
    {
        std::cout << "Cell index grid: " << std::endl;
        auto g = ccm->exterior_grid.cell_indices();
        auto s = g.shape();
        for(int i = 0; i < s[0]; ++i) {
            for(int j = 0; j < s[1]; ++j) {
                std::cout << g(i,j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << "Cell coordinates: " << ccm->exterior_grid.num_cells()  << std::endl;
    std::cout << "Cell coordinates: " << ccm->exterior_grid.cell_coords().size()  << std::endl;
    for(auto&& [a,b]: ccm->exterior_grid.cell_coords()) {
        std::cout << a << ":" << b << " ";
    }
    std::cout << std::endl;


    //set_colors(ccm.active_grid_cell_mask);

    {
        auto V = ccm->vertices();
        auto E = ccm->cut_edges_eigen();
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
        cutcell_mesh.setEdgeBuffer(V.cast<float>().eval(),E.cast<unsigned int>().eval());
        cutcell_drawable->deactivate();

        ccm->faces();
        grid_mesh.set(ccm->vertex_grid());
        grid_drawable->activate_edges();
    }
    update_faces();
}


void MeshViewer::update_face(int idx) {
    if(!ccm) return;

    if(idx < 0 || idx >= ccm->num_cells()) return;
    auto c = ccm->cell(idx);
    for(auto&& c: c) {
        std::copy
            (c.begin(),c.end(),std::ostream_iterator<int>(std::cout,","));
        std::cout << " ";
    }
    std::cout << std::endl;
    auto d = *c.begin();
    if(d.size() < 3) return;
    auto F = mtao::geometry::mesh::earclipping(ccm->vertices(),d);
    cutcell_face_mesh.setTriangleBuffer(F.cast<unsigned int>());
    cutcell_drawable->activate_triangles();
}

void MeshViewer::update_faces() {
    if(!ccm) return;

    if(colors.cols() != ccm->num_cells()) {
        colors.resize(4,ccm->num_cells());
        colors.setRandom();
        colors.row(3).setConstant(1);
    }
    std::vector<std::tuple<mtao::ColVecs3i, mtao::Vec4d>> FCs;
    for(int i = 0; i < ccm->num_cells(); ++i) {
        FCs.reserve(ccm->num_cells());
        auto c = ccm->cell(i);
        auto&& d = *c.begin();
        auto F = mtao::geometry::mesh::earclipping(ccm->vertices(),d);
        FCs.emplace_back(std::move(F), colors.col(i));
    }
    auto [V,F,C] = mtao::geometry::mesh::stack_meshes(ccm->vertices(), FCs);

    cutcell_mesh.setTriangleBuffer(V.cast<float>(),F.cast<unsigned int>());
    cutcell_mesh.setColorBuffer(C.cast<float>().eval());
    cutcell_drawable->activate_triangles();

}

MAGNUM_APPLICATION_MAIN(MeshViewer)

