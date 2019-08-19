ADD_LIBRARY(cutmesh_maker
    make_cutmesh_from_cmdline.cpp
    make_cutmesh_from_cmdline.hpp
    make_cutmesh_generator_from_cmdline.hpp
    )

TARGET_LINK_LIBRARIES(cutmesh_maker mandoline)
ADD_DEFINITIONS(-g)

ADD_EXECUTABLE(cutmesh_to_obj cutmesh_to_obj.cpp)
TARGET_LINK_LIBRARIES(cutmesh_to_obj cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao::common)
ADD_EXECUTABLE(write_to_obj write_to_obj.cpp)
TARGET_LINK_LIBRARIES(write_to_obj cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao::common)
ADD_EXECUTABLE(obj_to_cutmesh obj_to_cutmesh.cpp)
TARGET_LINK_LIBRARIES(obj_to_cutmesh cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao::common)


find_package(Corrade REQUIRED Utility)

ADD_EXECUTABLE(cfg_to_cutmesh cfg_to_cutmesh.cpp)
TARGET_LINK_LIBRARIES(cfg_to_cutmesh cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao::common  Corrade::Utility)

ADD_EXECUTABLE(remove_self_intersections remove_self_intersections.cpp
    ../include/mandoline/construction/remesh_self_intersections.hpp
    ../src/construction/remesh_self_intersections.cpp
    )
TARGET_LINK_LIBRARIES(remove_self_intersections gmp CGAL mpfr Boost::thread
    mtao::common igl::core)


ADD_EXECUTABLE(cutmesh_validation
    validation/cutmesh_validation.cpp 
    validation/faces_fully_utilized.cpp validation/paired_boundary.cpp validation/pcwn_count.cpp validation/possible_cells.cpp
    validation/grid_cells.cpp
    validation/region_volumes.cpp
    )

TARGET_LINK_LIBRARIES(cutmesh_validation gmp CGAL mpfr Boost::thread
    mandoline
    )

IF(${USE_OPENGL})
find_package(Magnum REQUIRED GL Shaders GlxContext GlfwApplication SceneGraph MeshTools)
find_package(MagnumIntegration REQUIRED ImGui)
ADD_EXECUTABLE(make_cutmesh_gui make_cutmesh_gui.cpp)
TARGET_LINK_LIBRARIES(make_cutmesh_gui mandoline OpenMP::OpenMP_CXX mtao::common mtao::visualization)

ENDIF(${USE_OPENGL})
