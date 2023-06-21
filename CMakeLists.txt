CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(RT CXX)

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF()

SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

SET(RT_SOURCES
        src/main.cpp )

SET(RT_INCLUDES
        include/vector.h
        include/ray.hpp
        include/hit.h
        include/plane.hpp
        include/triangle.hpp
        )

SET(CMAKE_CXX_STANDARD 14)



ADD_EXECUTABLE(${PROJECT_NAME} ${RT_SOURCES} ${RT_INCLUDES})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE include)
