CMAKE_MINIMUM_REQUIRED(VERSION 3.10.2)

project(oneMD DESCRIPTION "Molecular Dynamics simulator for Intel oneAPI." LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(DEFINED ENV{ONEAPI_ROOT})
  set (USE_ONEAPI ON)
  message(STATUS "Using Intel oneAPI environment at $ENV{ONEAPI_ROOT}.")
else()
  set (USE_ONEAPI OFF)
  message(WARNING "Not using Intel oneAPI. Only standard C++ features and libraries will be available.")
endif()

find_package (Threads REQUIRED)

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

set (XDRFILE_SRC
  ext/libxdrfile/src/trr_seek.c
  ext/libxdrfile/src/xtc_seek.c
  ext/libxdrfile/src/xdrfile_trr.c
  ext/libxdrfile/src/xdrfile_xtc.c
  ext/libxdrfile/src/xdrfile.c
)

set (CPP_SIMS_SRC src/simulators/JB.cpp)
if(USE_ONEAPI)
  set (SIMS_SRC ${CPP_SIMS_SRC})
else()
  set (SIMS_SRC ${CPP_SIMS_SRC})
endif()
set (SRC src/Simulator.cpp ${SIMS_SRC} src/main.cpp)

add_library (xdrfile SHARED ${XDRFILE_SRC} )
target_include_directories(xdrfile PRIVATE  
    ext/libxdrfile/include
)

add_executable(${TGT} ${EXT_SRC} ${SRC})
target_include_directories(${TGT} PRIVATE 
    ext/embedFiglet/src 
    ext/spdlog/include 
    ext/tclap-1.2.3/include 
    ext/libxdrfile/include
    src
)
target_link_libraries(${TGT} PUBLIC 
    ${CMAKE_THREAD_LIBS_INIT}
    xdrfile
)
target_compile_definitions(${TGT} PRIVATE SPDLOG_COMPILED_LIB=1)
target_compile_features(${TGT} PRIVATE cxx_std_17)

add_custom_command(
    TARGET ${TGT} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/${BIN}        
            ${CMAKE_SOURCE_DIR}/${BIN_OUT}
)

add_executable(jb_test ext/test_sims/jb_test.cpp)

MESSAGE( STATUS "C++ compiler                 = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "MPI                          = $ENV{I_MPI_ROOT}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE      = ${CMAKE_CXX_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_DEBUG        = ${CMAKE_CXX_FLAGS_DEBUG}" )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME            = ${CMAKE_SYSTEM_NAME}" )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR       = ${CMAKE_SYSTEM_PROCESSOR}" )
IF(APPLE)
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT          = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES    = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET= ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()
