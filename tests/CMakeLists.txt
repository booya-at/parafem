INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/headers)
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/vtkWriter/headers)

add_executable(truss1 truss1.cpp)
target_link_libraries(truss1 parafem vtkWriter)

add_executable(truss2 truss2.cpp)
target_link_libraries(truss2 parafem vtkWriter)

add_executable(membrane3_1 membrane3_1.cpp)
target_link_libraries(membrane3_1 parafem vtkWriter)

add_executable(membrane3_2 membrane3_2.cpp)
target_link_libraries(membrane3_2 parafem vtkWriter)

# normal spannungen
add_executable(membrane4_1 membrane4_1.cpp)
target_link_libraries(membrane4_1 parafem vtkWriter)

add_executable(membrane4_2 membrane4_2.cpp)
target_link_libraries(membrane4_2 parafem vtkWriter)

# scherspannungen
add_executable(membrane4_3 membrane4_3.cpp)
target_link_libraries(membrane4_3 parafem vtkWriter)


set(TESTS
	truss1
	truss2
	membrane3_1
	membrane3_2
	membrane4_1
	membrane4_2
	membrane4_3
)

install(TARGETS ${TESTS} DESTINATION share/parafem/tests)