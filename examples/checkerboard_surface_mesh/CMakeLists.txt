
ADD_EXECUTABLE(checkerboard_surface_mesh main.cpp)
TARGET_LINK_LIBRARIES(checkerboard_surface_mesh mandoline_cutmesh3 cxxopts
    spdlog)
