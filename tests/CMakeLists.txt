ADD_EXECUTABLE(adaptive_grid_test adaptive_grid_test.cpp)
TARGET_LINK_LIBRARIES(adaptive_grid_test mandoline_cutmesh3)

ADD_EXECUTABLE(projection_unity_test projection_unity_test.cpp)
TARGET_LINK_LIBRARIES(projection_unity_test mandoline_cutmesh3)

ADD_EXECUTABLE(cutface_axial_test cutface_axial_test.cpp)
TARGET_LINK_LIBRARIES(cutface_axial_test mandoline_cutmesh3)

ADD_EXECUTABLE(laplacian_test laplacian_test.cpp)
TARGET_LINK_LIBRARIES(laplacian_test mandoline_cutmesh3 OpenMP::OpenMP_CXX)


if(HANDLE_SELF_INTERSECTIONS)
FIND_PACKAGE(CGAL REQUIRED COMPONENTS Core)
ADD_EXECUTABLE(igl_boolean_performance_test igl_boolean_performance_test.cpp)
TARGET_LINK_LIBRARIES(igl_boolean_performance_test mtao::headers igl::core gmp
    CGAL::CGAL CGAL::CGAL_Core mtao::common)
    ENDIF()
ADD_EXECUTABLE(mandoline_construction_performance_test mandoline_construction_performance_test.cpp)
TARGET_LINK_LIBRARIES(mandoline_construction_performance_test cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao::common)

ADD_EXECUTABLE(single_cube_test single_cube_test.cpp)
TARGET_LINK_LIBRARIES(single_cube_test mandoline)
