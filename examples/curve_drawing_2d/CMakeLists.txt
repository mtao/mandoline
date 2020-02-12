ADD_EXECUTABLE(curve_drawing_2d
    plcurve2.cpp
    plcurve2.hpp
    #viewer2d.cpp
    #viewer2d.h
    viewer2d_main.cpp
    )
TARGET_LINK_LIBRARIES(curve_drawing_2d mandoline mtao::visualization)
