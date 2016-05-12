cmake_minimum_required(VERSION 2.6)
ADD_DEFINITIONS(-std=c++11)

project(parafem)

add_subdirectory("src/")
add_subdirectory("tests")