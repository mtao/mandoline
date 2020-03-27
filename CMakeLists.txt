CMAKE_MINIMUM_REQUIRED(VERSION 3.11.1)
project(Mandoline
    DESCRIPTION "A robust cut-cell mesh generator"
    LANGUAGES CXX)
SET(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(BUILD_TESTING AND CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

SET(MTAO_PATH "" CACHE STRING "location of mtao core libs" )
SET(LIBIGL_PATH "" CACHE STRING "location of libigl libs" )

OPTION(USE_OPENGL "Build opengl stuff" ON)
IF(APPLE)
    SET(USE_OPENMP_PRE OFF)
ELSE(APPLE)
    SET(USE_OPENMP_PRE ON)
ENDIF(APPLE)
OPTION(USE_OPENMP "Use OpenMP stuff" ${USE_OPENMP_PRE})
OPTION(HANDLE_SELF_INTERSECTIONS "Handle Self-Intersections" ON)
SET(MANDOLINE_BUILD_TESTING ${BUILD_TESTING})


include(mandoline_fetch_dependencies)


find_package(Eigen3 REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)



protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS proto/cutmesh.proto)
protobuf_generate_python(PROTO_PY proto/cutmesh.proto)

ADD_LIBRARY(cutmesh_proto OBJECT ${PROTO_SRCS} ${PROTO_HDRS})
TARGET_LINK_LIBRARIES(cutmesh_proto protobuf::libprotobuf)
TARGET_INCLUDE_DIRECTORIES(cutmesh_proto PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

IF(HANDLE_SELF_INTERSECTIONS)
    SET(MANDOLINE_SELF_INTERSECTION_LIBS
        gmp CGAL::CGAL mpfr Boost::thread
        )
ENDIF(HANDLE_SELF_INTERSECTIONS)

add_library(mandoline_headers INTERFACE)
target_include_directories(mandoline_headers INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include> ${CMAKE_CURRENT_BINARY_DIR}
)
TARGET_LINK_LIBRARIES(mandoline_headers PUBLIC INTERFACE Eigen3::Eigen
    mtao::headers mtao::geometry igl::core protobuf::libprotobuf
    )



SET(COMMON_SRCS
    #src/diffgeo_utils.cpp
    src/vertex.cpp
    src/domain_boundary.cpp
    src/operators/boundary.cpp
    src/interpolated_edge.cpp
    )


SET(CUTMESH2_SRCS
    ${COMMON_SRCS}
    src/mesh2.cpp
    src/operators/boundary2.cpp
    src/operators/interpolation2.cpp
    src/operators/volume2.cpp
    src/operators/diffgeo2.cpp
    src/cutface2.cpp
    )
SET(CUTMESH3_SRCS
    ${COMMON_SRCS}
    src/mesh3.cpp
    src/cutcell.cpp
    src/proto_util.cpp
    src/barycentric_triangle_face.cpp
    src/cutface3.cpp
    src/adaptive_grid.cpp
    src/operators/boundary3.cpp
    src/operators/diffgeo3.cpp
    src/operators/interpolation3.cpp
    src/operators/masks.cpp
    src/operators/volume3.cpp
    )



SET(CUTMESH_HDRS
    include/mandoline/line_impl.hpp
    include/mandoline/line.hpp
    include/mandoline/mesh_impl.hpp
    include/mandoline/mesh.hpp
    include/mandoline/coord_masked_geometry.hpp
    #include/mandoline/diffgeo_utils.hpp
    include/mandoline/coord_mask.hpp
    include/mandoline/adaptive_grid.hpp
    include/mandoline/exterior_grid.hpp
    include/mandoline/exterior_grid_impl.hpp
    include/mandoline/cutedge.hpp
    include/mandoline/cutedge_impl.hpp
    $<TARGET_OBJECTS:cutmesh_proto>
    )

SET(CUTMESH2_HDRS
    ${CUTMESH_HDRS}
    include/mandoline/mesh2.hpp
    include/mandoline/operators/boundary2.hpp
    )
SET(CUTMESH3_HDRS
    ${CUTMESH_HDRS}
    include/mandoline/proto_util.hpp
    include/mandoline/cutface.hpp
    include/mandoline/cutface_impl.hpp
    include/mandoline/barycentric_triangle_face.hpp
    include/mandoline/mesh3.hpp
    include/mandoline/cutcell.hpp
    include/mandoline/operators/boundary3.hpp
    include/mandoline/operators/interpolation3.hpp
    include/mandoline/operators/masks.hpp
    include/mandoline/operators/volume3.hpp
    )

SET(CONSTRUCTION_SRCS
    ${COMMON_SRCS}
    src/construction/generator3_make_faces.cpp
    src/construction/generator.cpp
    src/construction/generator3.cpp
    src/construction/cell_collapser.cpp
    src/construction/remesh_self_intersections.cpp
    src/construction/preprocess_mesh.cpp
    src/construction/generator2.cpp
    src/construction/subgrid_transformer.cpp
    src/construction/face_collapser.cpp
    src/construction/construct.cpp
    src/construction/adaptive_grid_factory.cpp
    )


SET(CONSTRUCTION_HDRS
    ${CUTMESH2_HDRS}
    ${CUTMESH3_HDRS}
    include/mandoline/construction/generator.hpp
    include/mandoline/construction/generator_impl.hpp
    include/mandoline/construction/vertex_types.hpp
    include/mandoline/construction/vertex_types_impl.hpp
    include/mandoline/construction/cutdata.hpp
    include/mandoline/construction/cutdata_impl.hpp
    include/mandoline/construction/facet_intersections.hpp
    include/mandoline/construction/facet_intersections_impl.hpp
    include/mandoline/construction/subgrid_transformer.hpp
    include/mandoline/construction/construct.hpp
    include/mandoline/construction/cell_collapser.hpp
    include/mandoline/construction/face_collapser.hpp
    include/mandoline/construction/adaptive_grid_factory.hpp
    include/mandoline/construction/remesh_self_intersections.hpp
    )

if(USE_OPENGL)
    set(CONSTRUCTION_SRCS ${CONSTRUCTION_SRCS} src/construction/construct_imgui.cpp)
    set(CONSTRUCTION_HDRS ${CONSTRUCTION_HDRS} src/construction/construct_imgui.cpp)
endif()

SET(MANDOLINE2_SRCS
    ${CUTMESH2_SRCS}
    ${CUTMESH2_HDRS}
    ${CONSTRUCTION_SRCS}
    ${CONSTRUCTION_HDRS}
    )
SET(MANDOLINE3_SRCS
    ${CUTMESH3_SRCS}
    ${CUTMESH3_HDRS}
    ${CONSTRUCTION_SRCS}
    ${CONSTRUCTION_HDRS}
    )

SET(MANDOLINE_CUTMESH_SRCS
    ${CUTMESH3_SRCS}
    ${CUTMESH3_HDRS}
    ${CONSTRUCTION_SRCS}
    ${CONSTRUCTION_HDRS}
    )


SET(TOOLS_SRCS
    src/tools/planar_slicer.cpp
    src/tools/exploded_mesh.cpp
    src/tools/cutmesh_info.cpp
    )
SET(TOOLS_HDRS
    include/mandoline/tools/planar_slicer.hpp
    include/mandoline/tools/exploded_mesh.hpp
    include/mandoline/tools/cutmesh_info.hpp
    )

SET(TOOLS2_SRCS
    src/tools/plcurve_io.cpp
    src/tools/edges_to_plcurves.cpp
    )
SET(TOOLS2_HDRS
    include/mandoline/tools/plcurve_io.hpp
    include/mandoline/tools/edges_to_plcurves.hpp
    )

ADD_LIBRARY(cutmesh2_obj OBJECT ${CUTMESH2_SRCS} ${CUTMESH2_HDRS})
ADD_LIBRARY(cutmesh3_obj OBJECT ${CUTMESH3_SRCS} ${CUTMESH3_HDRS})
ADD_LIBRARY(construction_obj OBJECT ${CONSTRUCTION_SRCS} ${CONSTRUCTION_HDRS})
ADD_LIBRARY(tools_obj OBJECT ${TOOLS_SRCS} ${TOOLS_HDRS})
ADD_LIBRARY(tools2_obj OBJECT ${TOOLS2_SRCS} ${TOOLS2_HDRS})



# TODO: figure out why transitivity doesnt work
TARGET_LINK_LIBRARIES(cutmesh2_obj PUBLIC mtao::headers mtao::geometry igl::core mandoline_headers protobuf::libprotobuf)
TARGET_LINK_LIBRARIES(cutmesh3_obj PUBLIC mtao::headers mtao::geometry igl::core mandoline_headers protobuf::libprotobuf)
TARGET_LINK_LIBRARIES(construction_obj PUBLIC mtao::headers mtao::geometry igl::core mandoline_headers protobuf::libprotobuf)

TARGET_LINK_LIBRARIES(tools_obj PUBLIC mtao::headers mtao::geometry igl::core mandoline_headers protobuf::libprotobuf)
TARGET_LINK_LIBRARIES(tools2_obj PUBLIC mandoline_headers mtao::headers mtao::geometry igl::core mandoline_headers protobuf::libprotobuf)

ADD_DEFINITIONS(-g)
IF(OpenMP_Found)
    ADD_DEFINITIONS(-DMTAO_OPENMP)
ENDIF(OpenMP_Found)


ADD_DEFINITIONS(-DUSE_FLOP_FREE_ANGLE_COMPUTATION)
ADD_LIBRARY(mandoline
    $<TARGET_OBJECTS:cutmesh_proto>
    $<TARGET_OBJECTS:cutmesh2_obj>
    $<TARGET_OBJECTS:cutmesh3_obj>
    $<TARGET_OBJECTS:construction_obj>
    $<TARGET_OBJECTS:tools_obj>
    $<TARGET_OBJECTS:tools2_obj>
    )
TARGET_LINK_LIBRARIES(mandoline
    mtao::geometry OpenMP::OpenMP_CXX
    ${MANDOLINE_SELF_INTERSECTION_LIBS}
    igl::core
    protobuf::libprotobuf
    mandoline_headers
    )

ADD_LIBRARY(mandoline2 EXCLUDE_FROM_ALL
    $<TARGET_OBJECTS:cutmesh_proto>
    $<TARGET_OBJECTS:cutmesh2_obj>
    $<TARGET_OBJECTS:construction_obj>
    $<TARGET_OBJECTS:tools2_obj>
    )
TARGET_LINK_LIBRARIES(mandoline2
    mtao::geometry OpenMP::OpenMP_CXX
    ${MANDOLINE_SELF_INTERSECTION_LIBS}
    igl::core
    protobuf::libprotobuf
    mandoline_headers
    )


ADD_LIBRARY(mandoline3 EXCLUDE_FROM_ALL
    $<TARGET_OBJECTS:cutmesh_proto>
    $<TARGET_OBJECTS:cutmesh3_obj>
    $<TARGET_OBJECTS:construction_obj>
    $<TARGET_OBJECTS:tools_obj>
    )
TARGET_LINK_LIBRARIES(mandoline3 
    mtao::geometry OpenMP::OpenMP_CXX
    ${MANDOLINE_SELF_INTERSECTION_LIBS} 
    igl::core
protobuf::libprotobuf
    mandoline_headers
    )

ADD_LIBRARY(mandoline_cutmesh3
    $<TARGET_OBJECTS:cutmesh_proto>
    $<TARGET_OBJECTS:cutmesh3_obj>
    $<TARGET_OBJECTS:tools_obj>
    )
TARGET_LINK_LIBRARIES(mandoline_cutmesh3
    mtao::geometry OpenMP::OpenMP_CXX
    ${MANDOLINE_SELF_INTERSECTION_LIBS}
    igl::core
    protobuf::libprotobuf
    mandoline_headers
    )
message(STATUS "BINARY DIR: ${CMAKE_CURRENT_BINARY_DIR}")

TARGET_INCLUDE_DIRECTORIES(mandoline PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
TARGET_INCLUDE_DIRECTORIES(mandoline2 PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
TARGET_INCLUDE_DIRECTORIES(mandoline3 PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
TARGET_INCLUDE_DIRECTORIES(mandoline_cutmesh3 PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

if(USE_OPENGL)
    TARGET_LINK_LIBRARIES(construction_obj PUBLIC ImGui::ImGui MagnumIntegration::ImGui)
endif()

IF(HANDLE_SELF_INTERSECTIONS)
    FOREACH(TARG mandoline mandoline2 mandoline3 construction_obj)
        TARGET_COMPILE_DEFINITIONS(${TARG}
            PUBLIC -DMANDOLINE_HANDLE_SELF_INTERSECTIONS)
    ENDFOREACH()
ENDIF(HANDLE_SELF_INTERSECTIONS)

ADD_SUBDIRECTORY(examples EXCLUDE_FROM_ALL)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND MANDOLINE_BUILD_TESTING)
    ADD_SUBDIRECTORY(tests)
endif()
ADD_SUBDIRECTORY(tools EXCLUDE_FROM_ALL)
