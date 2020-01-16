ADD_EXECUTABLE(slice_viewer main.cpp)
ADD_DEFINITIONS(-g)
TARGET_LINK_LIBRARIES(slice_viewer mandoline mtao::visualization)
