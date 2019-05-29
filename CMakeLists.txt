CMAKE_MINIMUM_REQUIRED(VERSION 3.10.1)
SET(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

ADD_SUBDIRECTORY(extern/core)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/core/extern/imgui
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/libigl/include
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/core/include)

FIND_PACKAGE(OpenMP)
FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(Boost COMPONENTS thread REQUIRED)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/unsupported)

find_package(Protobuf REQUIRED)
include_directories(${Protobuf_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS proto/cutmesh.proto)
protobuf_generate_python(PROTO_PY proto/cutmesh.proto)





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

ADD_LIBRARY(mandoline
    ${CUTMESH2_SRCS}
    ${CUTMESH3_HDRS}

    ${CUTMESH3_SRCS}
    ${CUTMESH2_HDRS}

    ${CONSTRUCTION_SRCS}
    ${CONSTRUCTION_HDRS}
    )
TARGET_LINK_LIBRARIES(mandoline
    mtao_geometry OpenMP::OpenMP_CXX
    gmp CGAL mpfr Boost::thread
    ${Protobuf_LIBRARIES}
    )

ADD_LIBRARY(mandoline2 EXCLUDE_FROM_ALL
    ${CUTMESH2_SRCS}
    ${CUTMESH2_HRDS}
    ${CONSTRUCTION_SRCS}
    ${CONSTRUCTION_HDRS}
    )
TARGET_LINK_LIBRARIES(mandoline2
    mtao_geometry OpenMP::OpenMP_CXX
    gmp CGAL mpfr Boost::thread
    ${Protobuf_LIBRARIES}
    )


ADD_LIBRARY(mandoline3 EXCLUDE_FROM_ALL
    ${CUTMESH3_SRCS}
    ${CUTMESH3_HRDS}
    ${CONSTRUCTION_SRCS}
    ${CONSTRUCTION_HDRS}
    )
TARGET_LINK_LIBRARIES(mandoline3 
    mtao_geometry OpenMP::OpenMP_CXX
    gmp CGAL mpfr Boost::thread
    ${Protobuf_LIBRARIES}
    )

ADD_LIBRARY(mandoline_cutmesh3
    ${CUTMESH3_SRCS}
    ${CUTMESH3_HRDS}
    )
TARGET_LINK_LIBRARIES(mandoline_cutmesh3
    mtao_geometry OpenMP::OpenMP_CXX
    gmp CGAL mpfr Boost::thread
    ${Protobuf_LIBRARIES}
    )


ADD_SUBDIRECTORY(examples)
