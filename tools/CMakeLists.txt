ADD_LIBRARY(cutmesh_maker
    make_cutmesh_from_cmdline.cpp
    make_cutmesh_from_cmdline.hpp
    make_cutmesh_generator_from_cmdline.hpp
    )

TARGET_LINK_LIBRARIES(cutmesh_maker mandoline cxxopts)
ADD_DEFINITIONS(-g)

ADD_EXECUTABLE(cutmesh_to_obj cutmesh_to_obj.cpp)
TARGET_LINK_LIBRARIES(cutmesh_to_obj cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao::common fmt)

ADD_EXECUTABLE(obj_to_cutmesh obj_to_cutmesh.cpp)
TARGET_LINK_LIBRARIES(obj_to_cutmesh cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao::common)

ADD_EXECUTABLE(cutmesh_info cutmesh_info.cpp)
TARGET_LINK_LIBRARIES(cutmesh_info mandoline_cutmesh3 nlohmann_json)

ADD_EXECUTABLE(add_triangulation_to_cutmesh add_triangulation_to_cutmesh.cpp)
TARGET_LINK_LIBRARIES(add_triangulation_to_cutmesh mandoline_cutmesh3 cxxopts)

ADD_EXECUTABLE(boundary_curves_to_cutmesh2 boundary_curves_to_cutmesh2.cpp)
TARGET_LINK_LIBRARIES(boundary_curves_to_cutmesh2 mandoline OpenMP::OpenMP_CXX mtao::common cxxopts)





if(MANDOLINE_USE_JSON)
ADD_EXECUTABLE(json_to_cutmesh json_to_cutmesh.cpp)
TARGET_LINK_LIBRARIES(json_to_cutmesh mandoline)

ADD_EXECUTABLE(json_to_cutmesh2 json_to_cutmesh2.cpp)
TARGET_LINK_LIBRARIES(json_to_cutmesh2 mandoline)
endif()
IF(Corrade_FOUND)
    ADD_EXECUTABLE(cfg_to_cutmesh cfg_to_cutmesh.cpp)
    TARGET_LINK_LIBRARIES(cfg_to_cutmesh cutmesh_maker mandoline OpenMP::OpenMP_CXX mtao::common  Corrade::Utility)
ENDIF(Corrade_FOUND)


IF(HANDLE_SELF_INTERSECTIONS)
    ADD_EXECUTABLE(remove_self_intersections remove_self_intersections.cpp
    $<TARGET_OBJECTS:construction_tools_obj>
        )
    TARGET_LINK_LIBRARIES(remove_self_intersections mandoline
        ${MANDOLINE_SELF_INTERSECTION_LIBS})
    TARGET_COMPILE_DEFINITIONS(remove_self_intersections
        PUBLIC -DMANDOLINE_HANDLE_SELF_INTERSECTIONS)


    ADD_LIBRARY(cutmesh_validation_tools
        validation/faces_fully_utilized.cpp validation/paired_boundary.cpp validation/pcwn_count.cpp validation/possible_cells.cpp
        validation/grid_cells.cpp
        validation/region_volumes.cpp
        validation/external_cell_valences.cpp
        validation/cutmesh_validation.cpp
        )
    TARGET_LINK_LIBRARIES(cutmesh_validation_tools gmp CGAL mpfr  fmt
        mandoline)

    if(DEFINED ENV{CIRCLECI})
        target_compile_definitions(cutmesh_validation_tools
            -DMANDOLINE_ASSUME_CELLS_DONT_INTERSECT)
    endif()

    ADD_EXECUTABLE(cutmesh_validator
        cutmesh_validator.cpp
        )

    TARGET_LINK_LIBRARIES(cutmesh_validator mandoline cutmesh_validation_tools
        Corrade::Utility
        nlohmann_json::nlohmann_json)
ENDIF(HANDLE_SELF_INTERSECTIONS)
if(${MANDOLINE_USE_OPENGL})
    ADD_EXECUTABLE(make_cutmesh_gui make_cutmesh_gui.cpp)
    TARGET_LINK_LIBRARIES(make_cutmesh_gui mandoline OpenMP::OpenMP_CXX mtao::common mtao::visualization)
    if(HANDLE_SELF_INTERSECTIONS) 
        TARGET_COMPILE_DEFINITIONS(remove_self_intersections
            PUBLIC -DMANDOLINE_HANDLE_SELF_INTERSECTIONS)
    TARGET_LINK_LIBRARIES(make_cutmesh_gui cutmesh_validation_tools)
    endif()

endif()

