ADD_LIBRARY(cutmesh_maker
    make_cutmesh_from_cmdline.cpp
    make_cutmesh_from_cmdline.hpp
    make_cutmesh_generator_from_cmdline.hpp
    )

TARGET_LINK_LIBRARIES(cutmesh_maker mandoline)
ADD_DEFINITIONS(-g)

ADD_EXECUTABLE(cutmesh_to_obj cutmesh_to_obj.cpp)
TARGET_LINK_LIBRARIES(cutmesh_to_obj cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao_core)
ADD_EXECUTABLE(write_to_obj write_to_obj.cpp)
TARGET_LINK_LIBRARIES(write_to_obj cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao_core)
ADD_EXECUTABLE(write_to_cutmesh write_to_cutmesh.cpp)
TARGET_LINK_LIBRARIES(write_to_cutmesh cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao_core)

IF(${USE_OPENGL})
find_package(Magnum REQUIRED GL Shaders GlxContext GlfwApplication SceneGraph MeshTools)
find_package(MagnumIntegration REQUIRED ImGui)
ADD_EXECUTABLE(make_cutmesh_gui make_cutmesh_gui.cpp)
TARGET_LINK_LIBRARIES(make_cutmesh_gui mandoline OpenMP::OpenMP_CXX mtao_core mtao_opengl)

ENDIF(${USE_OPENGL})
