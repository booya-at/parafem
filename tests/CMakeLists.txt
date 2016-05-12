include_directories( ${CMAKE_SOURCE_DIR}/src/headers)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_library(vtkwriter SHARED vtkWriter.cpp)
target_link_libraries(vtkwriter ${VTK_LIBRARIES} paraFEM)


add_executable(truss1 truss1.cpp)
target_link_libraries(truss1 paraFEM vtkwriter)

add_executable(truss2 truss2.cpp)
target_link_libraries(truss2 paraFEM vtkwriter)

add_executable(membrane3_1 membrane3_1.cpp)
target_link_libraries(membrane3_1 paraFEM vtkwriter)


add_executable(membrane3_2 membrane3_2.cpp)
target_link_libraries(membrane3_2 paraFEM vtkwriter)


add_executable(membrane4_1 membrane4_1.cpp)
target_link_libraries(membrane4_1 paraFEM vtkwriter)


add_executable(membrane4_2 membrane4_2.cpp)
target_link_libraries(membrane4_2 paraFEM vtkwriter)