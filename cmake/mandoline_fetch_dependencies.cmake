include(FetchContent REQUIRED)


set(MTAO_COMMIT 
    e895b0441588ddfa6b092150f38b3f5b9146a777
    )
set(PROTOBUF_COMMIT v3.11.3)
set(LIBIGL_COMMIT v2.1.0)
set(CATCH_COMMIT v2.9.1)
set(CGAL_COMMIT releases/CGAL-5.0.1)

function(fetch_dep REPO_NAME GIT_REPO GIT_TAG ADD_SUBDIR)
    FetchContent_Declare(
        ${REPO_NAME}
        GIT_REPOSITORY ${GIT_REPO}
        #GIT_TAG f6b406427400ed7ddb56cfc2577b6af571827c8c
        GIT_TAG ${GIT_TAG}
        )
    if(ADD_SUBDIR)
        if(${CMAKE_VERSION} VERSION_LESS 3.14)
            FetchContent_Populate(${REPO_NAME})
            add_subdirectory(${${REPO_NAME}_SOURCE_DIR} ${${REPO_NAME}_BINARY_DIR})
        else()
            FetchContent_MakeAvailable(${REPO_NAME})
        endif()
    else()
        FetchContent_Populate(${REPO_NAME})
    endif()
    set(${REPO_NAME}_SOURCE_DIR ${${REPO_NAME}_SOURCE_DIR} PARENT_SCOPE)
    set(${REPO_NAME}_BINARY_DIR ${${REPO_NAME}_BINARY_DIR} PARENT_SCOPE)
endfunction()

OPTION(MTAO_USE_ELTOPO "Should we build the el topo submodule" OFF)
OPTION(MTAO_USE_LOSTOPOS "Should we build the LOS Topos submodule" OFF)
OPTION(MTAO_USE_OPENGL "Build opengl stuff" ${USE_OPENGL})
OPTION(MTAO_USE_PNGPP "Use PNG++ for screenshots" OFF)

if(MTAO_PATH)
    ADD_SUBDIRECTORY("${MTAO_PATH}" ${CMAKE_BINARY_DIR}/mtao EXCLUDE_FROM_ALL)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${MTAO_PATH}/cmake")
else()

    fetch_dep(mtao https://github.com/mtao/core.git 
        ${MTAO_COMMIT}
        ON)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${mtao_SOURCE_DIR}/cmake")
endif()




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
    fetch_dep(libigl https://github.com/libigl/libigl.git ${LIBIGL_COMMIT} ON)
endif()

if(USE_OPENGL)
    find_package(ImGui REQUIRED)
    #find_package(MagnumIntegration COMPONENTS ImGui)
endif()

find_package(Protobuf)
if(NOT Protobuf_FOUND)
    option(protobuf_BUILD_TESTS "PROTOBUF BUILD TESTS" OFF)
    fetch_dep(protobuf https://github.com/protocolbuffers/protobuf.git ${PROTOBUF_COMMIT} ON)

endif()
if(HANDLE_SELF_INTERSECTIONS)
    find_package(CGAL COMPONENTS Core CONFIG REQUIRED)
    if(NOT CGAL_FOUND)
        find_package(CGAL COMPONENTS Core REQUIRED)
    endif()
    if(NOT CGAL_FOUND)
        #because its required we will never get here
        fetch_dep(
            CGAL
            https://github.com/CGAL/cgal.git
            ${CGAL_COMMIT}
            ON
            )

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

if(MANDOLINE_BUILD_TESTING)
    if(NOT Catch2_FOUND)
        if(NOT TARGET Catch2::Catch2)
        fetch_dep(
            catch2
            https://github.com/catchorg/Catch2.git
            ${CATCH_COMMIT}
            ON
            )
    endif()

    endif()
endif()
