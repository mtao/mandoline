CMAKE_MINIMUM_REQUIRED(VERSION 3.10.1)
SET(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


SET(MTAO_PATH "${CMAKE_CURRENT_SOURCE_DIR}/extern/core" CACHE STRING "location of mtao core libs" )
SET(LIBIGL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/extern/libigl" CACHE STRING "location of libigl libs" )
OPTION(USE_OPENGL "Build opengl stuff" ON)
OPTION(USE_OPENMP "Use OpenMP stuff" ON)
OPTION(HANDLE_SELF_INTERSECTIONS "Handle Self-Intersections" ON)

MESSAGE(STATUS "MTAO Path: ${MTAO_PATH}")
MESSAGE(STATUS "LIBIGL Path: ${LIBIGL_PATH}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${MTAO_PATH}/cmake")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${LIBIGL_PATH}/cmake")
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include
#    ${CMAKE_CURRENT_SOURCE_DIR}/extern/core/extern/imgui
#    ${CMAKE_CURRENT_SOURCE_DIR}/extern/libigl/include
#    ${CMAKE_CURRENT_SOURCE_DIR}/extern/core/include)


OPTION(MTAO_USE_ELTOPO "Should we build the el topo submodule" OFF)
OPTION(MTAO_USE_OPENGL "Build opengl stuff" ${USE_OPENGL})
OPTION(MTAO_USE_PNGPP "Use PNG++ for screenshots" OFF)
ADD_SUBDIRECTORY("${MTAO_PATH}" ${CMAKE_BINARY_DIR}/mtao EXCLUDE_FROM_ALL)


IF(USE_OPENMP)
    FIND_PACKAGE(OpenMP REQUIRED)
ENDIF(USE_OPENMP)

MESSAGE(STATUS "MODULE PATH:${CMAKE_MODULE_PATH}")
FIND_PACKAGE(Eigen3 REQUIRED)
#FIND_PACKAGE(libigl REQUIRED)
FIND_PACKAGE(Boost COMPONENTS thread REQUIRED)


option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         OFF)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PYTHOFF            "Use Python"                  OFF)
ADD_SUBDIRECTORY("${LIBIGL_PATH}" ${CMAKE_BINARY_DIR}/libigl EXCLUDE_FROM_ALL)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/unsupported)

find_package(Protobuf REQUIRED)
include_directories(${Protobuf_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS proto/cutmesh.proto)
protobuf_generate_python(PROTO_PY proto/cutmesh.proto)


find_package(CGAL)




SET(COMMON_SRCS
    #src/diffgeo_utils.cpp
    )


SET(CUTMESH2_SRCS
    src/mesh2.cpp
    )
SET(CUTMESH3_SRCS
    ${PROTO_SRCS} 
    src/mesh3.cpp
    src/cell.cpp
    src/proto_util.cpp
    src/barycentric_triangle_face.cpp
    src/cutface.cpp
    src/adaptive_grid.cpp
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
    )

SET(CUTMESH2_HDRS
    ${CUTMESH_HDRS}
    ${PROTO_HDRS}
    include/mandoline/mesh2.hpp
    )
SET(CUTMESH3_HDRS
    ${CUTMESH_HDRS}
    ${PROTO_HDRS}
    include/mandoline/construction/cell_collapser.hpp
    include/mandoline/construction/face_collapser.hpp
    include/mandoline/construction/remesh_self_intersections.hpp
    include/mandoline/proto_util.hpp
    include/mandoline/cutface.hpp
    include/mandoline/cutface_impl.hpp
    include/mandoline/barycentric_triangle_face.hpp
    include/mandoline/mesh3.hpp
    include/mandoline/cell.hpp
    )

SET(CONSTRUCTION_SRCS
    src/construction/generator3_make_faces.cpp
    src/construction/generator.cpp
    src/construction/generator3.cpp
    src/construction/cell_collapser.cpp
    src/construction/remesh_self_intersections.cpp
    src/construction/generator2.cpp
    src/construction/subgrid_transformer.cpp
    src/construction/face_collapser.cpp
    src/construction/construct.cpp
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
    )

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
    )
SET(TOOLS_HDRS
    include/mandoline/tools/planar_slicer.hpp
    include/mandoline/tools/exploded_mesh.hpp
    )

ADD_LIBRARY(cutmesh2_obj OBJECT ${CUTMESH2_SRCS} ${CUTMESH2_HDRS})
ADD_LIBRARY(cutmesh3_obj OBJECT ${CUTMESH3_SRCS} ${CUTMESH3_HDRS})
ADD_LIBRARY(construction_obj OBJECT ${CONSTRUCTION_SRCS} ${CONSTRUCTION_HDRS})
ADD_LIBRARY(tools_obj OBJECT ${TOOLS_SRCS} ${TOOLS_HDRS})

SET(MY_INCLUDES $<TARGET_PROPERTY:mtao::headers,INTERFACE_INCLUDE_DIRECTORIES>
 $<TARGET_PROPERTY:igl::core,INTERFACE_INCLUDE_DIRECTORIES>
 ${CMAKE_CURRENT_BINARY_DIR}
    )

MESSAGE(STATUS "CURRENT INCLUDES: ${MY_INCLUDES}")

TARGET_INCLUDE_DIRECTORIES(cutmesh2_obj PUBLIC ${MY_INCLUDES})
TARGET_INCLUDE_DIRECTORIES(cutmesh3_obj PUBLIC ${MY_INCLUDES})
TARGET_INCLUDE_DIRECTORIES(construction_obj PUBLIC ${MY_INCLUDES})
TARGET_INCLUDE_DIRECTORIES(tools_obj PUBLIC ${MY_INCLUDES})

ADD_DEFINITIONS(-g)
IF(OpenMP_Found)
    ADD_DEFINITIONS(-DMTAO_OPENMP)
ENDIF(OpenMP_Found)

IF(HANDLE_SELF_INTERSECTIONS)
    SET(SELF_INTERSECTION_LIB
        gmp CGAL::CGAL mpfr
        )
ENDIF(HANDLE_SELF_INTERSECTIONS)

ADD_DEFINITIONS(-DUSE_FLOP_FREE_ANGLE_COMPUTATION)
ADD_LIBRARY(mandoline
    $<TARGET_OBJECTS:cutmesh2_obj>
    $<TARGET_OBJECTS:cutmesh3_obj>
    $<TARGET_OBJECTS:construction_obj>
    $<TARGET_OBJECTS:tools_obj>
    )
TARGET_LINK_LIBRARIES(mandoline
    mtao::geometry OpenMP::OpenMP_CXX
    ${SELF_INTERSECTION_LIBS} Boost::thread
    igl::core
    ${Protobuf_LIBRARIES}
    )

ADD_LIBRARY(mandoline2 EXCLUDE_FROM_ALL
    $<TARGET_OBJECTS:cutmesh2_obj>
    $<TARGET_OBJECTS:construction_obj>
    )
TARGET_LINK_LIBRARIES(mandoline2
    mtao::geometry OpenMP::OpenMP_CXX
    ${SELF_INTERSECTION_LIBS} Boost::thread
    igl::core
    ${Protobuf_LIBRARIES}
    )


ADD_LIBRARY(mandoline3 EXCLUDE_FROM_ALL
    $<TARGET_OBJECTS:cutmesh3_obj>
    $<TARGET_OBJECTS:construction_obj>
    $<TARGET_OBJECTS:tools_obj>
    )
TARGET_LINK_LIBRARIES(mandoline3 
    mtao::geometry OpenMP::OpenMP_CXX
    ${SELF_INTERSECTION_LIBS} Boost::thread
    igl::core
    ${Protobuf_LIBRARIES}
    )

ADD_LIBRARY(mandoline_cutmesh3
    $<TARGET_OBJECTS:cutmesh3_obj>
    )
TARGET_LINK_LIBRARIES(mandoline_cutmesh3
    mtao::geometry OpenMP::OpenMP_CXX
    ${SELF_INTERSECTION_LIBS} Boost::thread
    igl::core
    ${Protobuf_LIBRARIES}
    )

IF(HANDLE_SELF_INTERSECTIONS)
    FOREACH(TARG mandoline mandoline2 mandoline3)
        TARGET_COMPILE_DEFINITIONS(${TARG}
            PUBLIC MANDOLINE_HANDLE_SELF_INTERSECTIONS)
    ENDFOREACH()
ENDIF(HANDLE_SELF_INTERSECTIONS)

ADD_SUBDIRECTORY(examples EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(tests EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(tools EXCLUDE_FROM_ALL)
