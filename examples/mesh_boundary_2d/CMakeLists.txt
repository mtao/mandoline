ADD_EXECUTABLE(mesh_boundary_2d
    main.cpp
    )
TARGET_LINK_LIBRARIES(mesh_boundary_2d mandoline mtao::visualization
    )
