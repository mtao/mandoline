#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/geometry/bounding_box.hpp>

//#include "mandoline/convex_gerrymanderer.hpp"
#include <mandoline/mesh3.hpp>
#include <mandoline/operators/boundary3.hpp>
#include <mandoline/operators/centroids3.hpp>
#include <mandoline/operators/volume3.hpp>

using namespace mandoline;

CutCellMesh<3> load_from_file(const std::string& filename) {
    return CutCellMesh<3>::from_file(filename);
}

PYBIND11_MODULE(mandoline_py, m) {
    m.doc() = "pybind11 example plugin";  // optional module docstring

//    pybind11::class_<CutCellMesh<3>>(m, "CutCellMesh")
//        .def("vertices", &CutCellMesh<3>::vertices)
//        .def("origV", &CutCellMesh<3>::origV, "vertices from the input mesh")
//        .def("origE", &CutCellMesh<3>::origE, "edges from the input mesh")
//        .def("origF", &CutCellMesh<3>::origF, "faces from the input mesh")
//        .def("triangulate_face", &CutCellMesh<3>::triangulate_face)
//        .def("subVs", &CutCellMesh<3>::compute_subVs, "Projects vertices to each axis-aligned plane; useful for doing 2D operations on axis-aligned cut-faces")
//        .def("triangulated_cell", &CutCellMesh<3>::triangulated_cell, "A triangulation of a particular cell")
//        .def("triangulate_faces", &CutCellMesh<3>::triangulate_faces, "Caches a triangulation of the faces of the mesh")
//        .def("regions", &CutCellMesh<3>::regions);
//    ;
//    m.def("load_from_file", &load_from_file,
//          "Load a cut-cell mesh from a cutmesh file");
//    m.def("cell_volumes", &operators::cell_volumes,
//          "The volume of each cell");
//    m.def("face_volumes", &operators::face_volumes,
//          "The volume of each face");
//    m.def("boundary", &operators::boundary,
//          "Boundary matrix for each cell");
}
