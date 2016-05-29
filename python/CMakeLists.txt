include_directories( ${CMAKE_SOURCE_DIR}/src/headers)

set(Python_ADDITIONAL_VERSIONS ${py})
find_package(PythonInterp ${py})
find_package(PythonLibs ${py})

EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c 
  "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
  OUTPUT_VARIABLE python_libs OUTPUT_STRIP_TRAILING_WHITESPACE )


include_directories( ${PYTHON_INCLUDE_DIRS})    # python.h

PYTHON_ADD_MODULE(paraFEM_py SHARED paraFEM_py.cpp)
target_link_libraries(paraFEM_py paraFEM)