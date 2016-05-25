cmake_minimum_required(VERSION 2.6)
ADD_DEFINITIONS(-std=c++11)

project(parafem)


find_package(OpenMP)

if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
endif()

add_subdirectory("src/")
add_subdirectory("tests")