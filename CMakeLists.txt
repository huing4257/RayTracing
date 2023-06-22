CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(RT CXX)

SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb -fopenmp")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -fopenmp")

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF()

ADD_SUBDIRECTORY(deps/vecmath)

SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

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

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)


ADD_EXECUTABLE(${PROJECT_NAME} ${RT_SOURCES} ${RT_INCLUDES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} vecmath)
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE include)
