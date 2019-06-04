ADD_EXECUTABLE(slice_viewer main.cpp slice_generator.cpp)
find_package(Magnum REQUIRED GL Shaders GlxContext GlfwApplication SceneGraph MeshTools)
find_package(MagnumIntegration REQUIRED ImGui)
TARGET_LINK_LIBRARIES(slice_viewer mandoline mtao_opengl)
