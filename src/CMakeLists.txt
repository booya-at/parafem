include_directories(${CMAKE_SOURCE_DIR}/src/headers)

set(HEADER_FILES ${CMAKE_SOURCE_DIR}/src/headers/base.h
                 ${CMAKE_SOURCE_DIR}/src/headers/case.h
                 ${CMAKE_SOURCE_DIR}/src/headers/element.h
                 ${CMAKE_SOURCE_DIR}/src/headers/node.h
                 ${CMAKE_SOURCE_DIR}/src/headers/material.h
                 ${CMAKE_SOURCE_DIR}/src/headers/constraint.h)

set(CPP_FILES  ${CMAKE_SOURCE_DIR}/src/case.cpp
               ${CMAKE_SOURCE_DIR}/src/node.cpp
               ${CMAKE_SOURCE_DIR}/src/element.cpp
               ${CMAKE_SOURCE_DIR}/src/truss.cpp
               ${CMAKE_SOURCE_DIR}/src/membrane.cpp
               ${CMAKE_SOURCE_DIR}/src/material.cpp)

add_library(paraFEM STATIC ${CPP_FILES} ${HEADER_FILES})

install(TARGETS paraFEM DESTINATION lib)

if (BUILD_VTK_WRITER)
    ADD_DEFINITIONS(-DBUILD_VTK_WRITER)
    ADD_SUBDIRECTORY("vtkWriter/")
ENDIF(BUILD_VTK_WRITER)

IF (BUILD_PYTHON_BINDINGS)
    ADD_SUBDIRECTORY("python/")
ENDIF(BUILD_PYTHON_BINDINGS)