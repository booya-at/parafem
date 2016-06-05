include_directories(${CMAKE_SOURCE_DIR}/src/headers)

set(HEADER_FILES ${CMAKE_SOURCE_DIR}/src/headers/base.h
                 ${CMAKE_SOURCE_DIR}/src/headers/case.h
                 ${CMAKE_SOURCE_DIR}/src/headers/element.h
                 ${CMAKE_SOURCE_DIR}/src/headers/node.h
                 ${CMAKE_SOURCE_DIR}/src/headers/material.h)

set(CPP_FILES  ${CMAKE_SOURCE_DIR}/src/case.cpp
               ${CMAKE_SOURCE_DIR}/src/node.cpp
               ${CMAKE_SOURCE_DIR}/src/element.cpp
               ${CMAKE_SOURCE_DIR}/src/truss.cpp
               ${CMAKE_SOURCE_DIR}/src/membrane.cpp
               ${CMAKE_SOURCE_DIR}/src/material.cpp)

add_library(paraFEM SHARED ${CPP_FILES} ${HEADER_FILES})

install(TARGETS paraFEM DESTINATION lib)

ADD_SUBDIRECTORY("vtkWriter/")

IF (BUILD_PYTHON_BINDINGS)
    ADD_SUBDIRECTORY("python/")
ENDIF(BUILD_PYTHON_BINDINGS)