
ADD_EXECUTABLE(particle_cell_ownership_benchmark main.cpp)
TARGET_LINK_LIBRARIES(particle_cell_ownership_benchmark mandoline
    mtao::common
    cxxopts)
