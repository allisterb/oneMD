CMAKE_MINIMUM_REQUIRED(VERSION 3.10.2)

project(oneMD DESCRIPTION "Molecular Dynamics simulator for Intel oneAPI." LANGUAGES CXX)

find_package (Threads REQUIRED)

add_subdirectory(ext/spdlog)

set (FIGLET 
    ext/embedFiglet/src/Figlet.hh 
    ext/embedFiglet/src/Figlet.cc
    ext/embedFiglet/src/Figlet_Font_standard.cc
    ext/embedFiglet/src/Figlet_Font_small.cc
    ext/spdlog/src/async.cpp
    ext/spdlog/src/cfg.cpp
    ext/spdlog/src/color_sinks.cpp
    ext/spdlog/src/file_sinks.cpp
    ext/spdlog/src/fmt.cpp
    ext/spdlog/src/spdlog.cpp
    ext/spdlog/src/stdout_sinks.cpp
)
set (SOURCE src/main.cpp)
add_executable(main ${FIGLET} ${SOURCE})
target_include_directories(main PRIVATE ext/embedFiglet/src ext/spdlog/include)
target_compile_definitions(main PRIVATE SPDLOG_COMPILED_LIB=1)
target_link_libraries(main  PUBLIC ${CMAKE_THREAD_LIBS_INIT})

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
