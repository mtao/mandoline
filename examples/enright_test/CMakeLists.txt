CMAKE_MINIMUM_REQUIRED(VERSION 3.7)
ADD_EXECUTABLE(enright_test
    src/main.cpp
    include/enright.h
    )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
TARGET_LINK_LIBRARIES(enright_test mtao::visualization mtao::geometry)

