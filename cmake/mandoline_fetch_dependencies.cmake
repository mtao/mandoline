include(FetchContent REQUIRED)

MESSAGE(STATUS "MANDOLINE FETCH DEPS")
MESSAGE(STATUS "MTAO Path: ${MTAO_PATH}")


OPTION(MTAO_USE_ELTOPO "Should we build the el topo submodule" OFF)
OPTION(MTAO_USE_LOSTOPOS "Should we build the LOS Topos submodule" OFF)
OPTION(MTAO_USE_OPENGL "Build opengl stuff" ${USE_OPENGL})
OPTION(MTAO_USE_PNGPP "Use PNG++ for screenshots" OFF)

if(MTAO_PATH)
    ADD_SUBDIRECTORY("${MTAO_PATH}" ${CMAKE_BINARY_DIR}/mtao EXCLUDE_FROM_ALL)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${MTAO_PATH}/cmake")
else()

    FetchContent_Declare(
        mtao
        GIT_REPOSITORY https://github.com/mtao/core.git
        GIT_TAG 0a57f1cf30e272f4699add60d035172a63acebd4
        )
    if(${CMAKE_VERSION} VERSION_LESS 3.14)
        FetchContent_Populate(mtao)
        add_subdirectory(${mtao_SOURCE_DIR} ${mtao_BINARY_DIR} EXCLUDE_FROM_ALL)
    else()
        FetchContent_MakeAvailable(mtao)
    endif()
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${mtao_SOURCE_DIR}/cmake")
endif()


FIND_PACKAGE(Eigen3 3.3.9)

if(NOT Eigen3_FOUND)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG bcbaad6d874d451817457ae0603f953cda3c0c06 
        )
    if(${CMAKE_VERSION} VERSION_LESS 3.14)
        FetchContent_Populate(eigen)
        add_subdirectory(${eigen_SOURCE_DIR} ${eigen_BINARY_DIR} EXCLUDE_FROM_ALL)
    else()
        FetchContent_MakeAvailable(eigen)
    endif()

ENDIF()

IF(USE_OPENMP)
    FIND_PACKAGE(OpenMP REQUIRED)
ENDIF(USE_OPENMP)

MESSAGE(STATUS "MODULE PATH:${CMAKE_MODULE_PATH}")
#FIND_PACKAGE(libigl REQUIRED)
FIND_PACKAGE(Boost COMPONENTS thread REQUIRED)


MESSAGE(STATUS "LIBIGL Path: ${LIBIGL_PATH}")
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
option(LIBIGL_SKIP_DOWNLOAD "Skip downloading external libraries" ON)
if(LIBIGL_PATH)
    ADD_SUBDIRECTORY("${LIBIGL_PATH}" ${CMAKE_BINARY_DIR}/libigl EXCLUDE_FROM_ALL)
else()
    FetchContent_Declare(
        libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        #GIT_TAG f6b406427400ed7ddb56cfc2577b6af571827c8c
        GIT_TAG v2.1.0
        )

    if(${CMAKE_VERSION} VERSION_LESS 3.14)
        FetchContent_Populate(libigl)
        add_subdirectory(${libigl_SOURCE_DIR} ${libigl_BINARY_DIR} EXCLUDE_FROM_ALL)
    else()
        FetchContent_MakeAvailable(libigl)
    endif()
endif()

if(USE_OPENGL)
    find_package(ImGui REQUIRED)
    find_package(MagnumIntegration COMPONENTS ImGui)
endif()

find_package(Protobuf REQUIRED)
find_package(CGAL)
