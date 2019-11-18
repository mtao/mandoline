ADD_EXECUTABLE(exploded_mesh main.cpp)
find_package(Magnum REQUIRED GL Shaders GlxContext GlfwApplication SceneGraph
    MeshTools Trade)
find_package(MagnumIntegration REQUIRED ImGui)
TARGET_LINK_LIBRARIES(exploded_mesh mandoline mtao::visualization)
