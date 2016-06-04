set(Python_ADDITIONAL_VERSIONS ${py})
find_package(PythonInterp ${py})
find_package(PythonLibs ${py})

EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c 
  "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
  OUTPUT_VARIABLE python_libs OUTPUT_STRIP_TRAILING_WHITESPACE )

SET(PYTHON_INSTALL_DIR ${python_libs}/paraFEM)

INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS})    # python.h, pybind11
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/headers)
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/vtkWriter/headers)

add_library(_paraFEM SHARED _paraFEM.cpp)
target_link_libraries(_paraFEM paraFEM vtkWriter)
set_target_properties(_paraFEM PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/paraFEM)
set_target_properties(_paraFEM PROPERTIES PREFIX "")

configure_file(test.py ${CMAKE_BINARY_DIR}/paraFEM/test.py COPYONLY)
configure_file(__init__.py ${CMAKE_BINARY_DIR}/paraFEM/__init__.py COPYONLY)


install(TARGETS _paraFEM DESTINATION ${PYTHON_INSTALL_DIR})
install(DIRECTORY "" DESTINATION ${PYTHON_INSTALL_DIR} FILES_MATCHING PATTERN "*.py")