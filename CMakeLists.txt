CMAKE_MINIMUM_REQUIRED(VERSION 3.10.2)

project(oneMD DESCRIPTION "Molecular Dynamics simulator for Intel oneAPI." LANGUAGES CXX)

find_package (Threads REQUIRED)

add_subdirectory(ext/spdlog)

if (CMAKE_BUILD_TYPE MATCHES Debug)
  set(TGT omdd)
  if (MSVC)
    set(BIN Debug/${TGT}.exe)
    set(BIN_OUT ${TGT}.exe)
  else()
    set (BIN ${TGT})
    set(BIN_OUT ${TGT})
  endif()
else ()
  set(TGT omd)
  if (MSVC)
    set(BIN Release/${TGT}.exe)
    set(BIN_OUT ${TGT}.exe)
  else()
    set (BIN ${TGT})
    set(BIN_OUT ${TGT})
  endif()
endif() 

set (EXT_SRC 
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
set (SRC src/main.cpp src/Simulator.hh)
add_executable(${TGT} ${EXT_SRC} ${SRC})
target_include_directories(${TGT} PRIVATE ext/embedFiglet/src ext/spdlog/include ext/tclap-1.4.0-rc1/include)
target_compile_definitions(${TGT} PRIVATE SPDLOG_COMPILED_LIB=1)
target_link_libraries(${TGT} PUBLIC ${CMAKE_THREAD_LIBS_INIT})

add_custom_command(
    TARGET ${TGT} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/${BIN}        
            ${CMAKE_SOURCE_DIR}/${BIN_OUT}
)
MESSAGE( STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_DEBUG       = ${CMAKE_CXX_FLAGS_DEBUG}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}" )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}" )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}" )
IF(APPLE)
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()
