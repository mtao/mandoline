include(FetchContent REQUIRED)



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
        GIT_TAG 08fa08710134a634ebbd6b4e7fea80a6a6ac4273
        )
    if(${CMAKE_VERSION} VERSION_LESS 3.14)
        FetchContent_Populate(mtao)
        add_subdirectory(${mtao_SOURCE_DIR} ${mtao_BINARY_DIR} EXCLUDE_FROM_ALL)
    else()
        FetchContent_MakeAvailable(mtao)
    endif()
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${mtao_SOURCE_DIR}/cmake")
endif()



if(NOT Eigen3_FOUND)
    FIND_PACKAGE(Eigen3 3.3.9)
ENDIF()
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
    #find_package(MagnumIntegration COMPONENTS ImGui)
endif()

find_package(Protobuf)
if(NOT Protobuf_FOUND)
    option(protobuf_BUILD_TESTS "PROTOBUF BUILD TESTS" OFF)
    FetchContent_Declare(
        protobuf
        GIT_REPOSITORY https://github.com/protocolbuffers/protobuf.git
        #GIT_TAG f6b406427400ed7ddb56cfc2577b6af571827c8c
        GIT_TAG v3.11.3
        )
        FetchContent_Populate(protobuf)
        add_subdirectory(${protobuf_SOURCE_DIR}/cmake ${protobuf_BINARY_DIR} EXCLUDE_FROM_ALL)

endif()
if(HANDLE_SELF_INTERSECTIONS)
    find_package(CGAL COMPONENTS Core CONFIG REQUIRED)
    if(NOT CGAL_FOUND)
        find_package(CGAL COMPONENTS Core REQUIRED)
    endif()
    if(NOT CGAL_FOUND)
        #because its required we will never get here
        FetchContent_Declare(
            CGAL
            GIT_REPOSITORY https://github.com/CGAL/cgal.git
            #GIT_TAG f6b406427400ed7ddb56cfc2577b6af571827c8c
            GIT_TAG releases/CGAL-5.0.1
            )

        if(${CMAKE_VERSION} VERSION_LESS 3.14)
            FetchContent_Populate(CGAL)
            add_subdirectory(${CGAL_SOURCE_DIR} ${CGAL_BINARY_DIR} EXCLUDE_FROM_ALL)
        else()
            FetchContent_MakeAvailable(CGAL)
        endif()
    endif(NOT CGAL_FOUND)
    FIND_PACKAGE(Boost COMPONENTS thread)
#if(NOT Boost_FOUND)
#    FetchContent_Declare(
#        boost
#        GIT_REPOSITORY https://github.com/boostorg/boost.git
#        GIT_TAG 789155d 
#        )
#    if(${CMAKE_VERSION} VERSION_LESS 3.14)
#        FetchContent_Populate(boost)
#        add_subdirectory(${boost_SOURCE_DIR} ${boost_BINARY_DIR} EXCLUDE_FROM_ALL)
#    else()
#        FetchContent_MakeAvailable(boost)
#    endif()
#     set(Boost_INCLUDE_DIR ${boost_INCLUDE_DIR} CACHE STRING "The location of hte boost include")
#     #set(Boost_NO_BOOST_CMAKE ON CACHE BOOL "We fetched our own boost so tell cmake to not look for it anymore")


#ENDIF()
endif(HANDLE_SELF_INTERSECTIONS)
