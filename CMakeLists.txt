CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(RT CXX)

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF()

find_package(OpenMP )
if(OpenMP_FOUND)
    message(STATUS "found openmp")
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS}, ${OPENMP_C_FLAGS})
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}, ${OPENMP_CXX_FLAGS})
endif()

SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

ADD_SUBDIRECTORY(deps/vecmath)

SET(RT_SOURCES
        src/main.cpp
        src/mesh.cpp
        )

SET(RT_INCLUDES
        include/ray.hpp
        include/hit.h
        include/plane.hpp
        include/triangle.hpp
        include/sphere.hpp
        include/mesh.hpp
        include/curve.hpp include/revsurface.hpp include/texture.h)

SET(CMAKE_CXX_STANDARD 14)



ADD_EXECUTABLE(${PROJECT_NAME} ${RT_SOURCES} ${RT_INCLUDES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} vecmath)
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE include)
