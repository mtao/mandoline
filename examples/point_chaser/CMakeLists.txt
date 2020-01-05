IF(USE_OPENGL)
ADD_EXECUTABLE(point_chaser main.cpp )
TARGET_LINK_LIBRARIES(point_chaser mandoline mtao::visualization)
ENDIF(USE_OPENGL)
