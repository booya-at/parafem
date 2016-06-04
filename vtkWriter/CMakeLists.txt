include_directories( ${CMAKE_SOURCE_DIR}/src/headers)
include_directories( ${CMAKE_SOURCE_DIR}/vtkWriter/headers)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_library(vtkWriter SHARED vtkWriter.cpp)
target_link_libraries(vtkWriter ${VTK_LIBRARIES} paraFEM)

install(TARGETS vtkWriter DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)