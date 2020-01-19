set(Python_ADDITIONAL_VERSIONS ${py})
find_package(PythonInterp ${py})
find_package(PythonLibs ${py})

find_package(pybind11 REQUIRED)

EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c 
  "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
  OUTPUT_VARIABLE python_libs OUTPUT_STRIP_TRAILING_WHITESPACE )

SET(PYTHON_INSTALL_DIR ${python_libs}/parafem)

INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS})    # python.h, pybind11
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/headers)

add_library(_parafem SHARED _parafem.cpp)
	target_link_libraries(_parafem parafem)

if (BUILD_VTK_WRITER)
	INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/vtkWriter/headers)
	target_link_libraries(_parafem vtkWriter)
endif (BUILD_VTK_WRITER)

set_target_properties(_parafem PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/parafem)
# !!!!!!!!!remember
set_target_properties(_parafem PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set_target_properties(_parafem PROPERTIES PREFIX "")
if( WIN32 )
    set_target_properties(_paraBEM PROPERTIES SUFFIX ".pyd")
endif()

configure_file(__init__.py ${CMAKE_BINARY_DIR}/parafem/__init__.py COPYONLY)


install(TARGETS _parafem DESTINATION ${PYTHON_INSTALL_DIR})
install(DIRECTORY "" DESTINATION ${PYTHON_INSTALL_DIR} FILES_MATCHING PATTERN "*.py")