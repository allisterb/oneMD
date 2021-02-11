CMAKE_MINIMUM_REQUIRED(VERSION 3.10.2)

project(oneMD VERSION 0.1 DESCRIPTION "Molecular Dynamics simulator for Intel oneAPI." LANGUAGES CXX)

add_subdirectory(ext/embedFiglet)

add_library(embedFiglet STATIC IMPORTED)
set_target_properties(embedFiglet PROPERTIES
    IMPORTED_LOCATION "ext/embedFiglet/Debug/embedFiglet_win__static.lib"
)

set(HEADERS ext/embedFiglet/src/Figlet.hh)
set (SOURCE src/main.cpp)
add_executable(main ${HEADERS} ${SOURCE})
target_include_directories(main PUBLIC ext/embedFiglet/src)
target_link_libraries(main embedFiglet)
MESSAGE( STATUS "Using ${SSE_FLAGS} extensions")
MESSAGE( STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}" )
MESSAGE( STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}" )
MESSAGE( STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}" )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}" )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}" )
IF(APPLE)
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()
MESSAGE( STATUS "BUILD_EXECUTABLE            = ${BUILD_EXECUTABLE}" )
