CMAKE_MINIMUM_REQUIRED(VERSION 3.9)
PROJECT(RT CXX)

ADD_SUBDIRECTORY(deps/vecmath)

# SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


SET(RT_SOURCES
        src/mesh.cpp
        src/main.cpp
        )

SET(RT_INCLUDES
        include/ray.hpp
        include/hit.h
        include/plane.hpp
        include/triangle.hpp
        include/sphere.hpp
        include/mesh.hpp
        include/curve.hpp include/revsurface.hpp include/texture.h include/moving.h)

SET(CMAKE_CXX_STANDARD 14)


ADD_EXECUTABLE(${PROJECT_NAME} ${RT_SOURCES} ${RT_INCLUDES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

TARGET_LINK_LIBRARIES(${PROJECT_NAME} PUBLIC vecmath)
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE include)
