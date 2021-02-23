CMAKE_MINIMUM_REQUIRED(VERSION 3.10.2)

project(oneMD DESCRIPTION "Data-parallel molecular dynamics simulator for Intel oneAPI." LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS 1)

if(DEFINED ENV{ONEAPI_ROOT})
  set(USE_ONEAPI ON)
else()
  set(USE_ONEAPI OFF)
endif()  

if (USE_ONEAPI)
  message(STATUS "Using Intel oneAPI environment at $ENV{ONEAPI_ROOT}.")
  if(WIN32)
    set(CMAKE_CXX_COMPILER "dpcpp-cl")
  else()
    set(CMAKE_CXX_COMPILER "dpcpp")
  endif()
else()
  message(WARNING "Not using Intel oneAPI. Only standard C++ features and libraries will be available.")
  find_package(OpenMP REQUIRED)
endif()

find_package (Threads REQUIRED)
find_package(Boost)

set (SPDLOG_SRC 
    ext/spdlog/src/async.cpp
    ext/spdlog/src/cfg.cpp
    ext/spdlog/src/color_sinks.cpp
    ext/spdlog/src/file_sinks.cpp
    ext/spdlog/src/fmt.cpp
    ext/spdlog/src/spdlog.cpp
    ext/spdlog/src/stdout_sinks.cpp
)
add_library (spdlog STATIC ${SPDLOG_SRC})
target_include_directories(spdlog PRIVATE ext/spdlog/include)
target_compile_definitions(spdlog PRIVATE SPDLOG_COMPILED_LIB=1)

set (FIGLET_SRC 
    ext/embedFiglet/src/Figlet.cc
    ext/embedFiglet/src/Figlet_Font_standard.cc
    ext/embedFiglet/src/Figlet_Font_small.cc
)

set (XDRFILE_SRC
  ext/libxdrfile/src/trr_seek.c
  ext/libxdrfile/src/xtc_seek.c
  ext/libxdrfile/src/xdrfile_trr.c
  ext/libxdrfile/src/xdrfile_xtc.c
  ext/libxdrfile/src/xdrfile.c
)
add_library (xdrfile STATIC ${XDRFILE_SRC} )
target_include_directories(xdrfile PRIVATE  
    ext/libxdrfile/include
)
if(MSVC)
	target_compile_definitions(xdrfile PRIVATE
		ftello=_ftelli64
		fseeko=_fseeki64
		off_t=long
	)
endif()

set (INIH_SRC ext/inih/ini.c)
add_library (inih STATIC ${INIH_SRC})
target_include_directories(inih PRIVATE ext/inih)

set (SIMS_SRC src/simulators/JB.cpp src/simulators/LJ.cpp)

set (SRC 
    ext/inih/cpp/INIReader.cpp
    src/Simulator.cpp 
    src/Vector.cpp 
    src/CubicBox.cpp 
    src/NeighborList.cpp 
    src/PdbFile.cpp
    src/Rdf.cpp 
    src/ThermodynamicVariable.cpp
    src/Thermostat.cpp
    src/Velocity.cpp
    ${SIMS_SRC} 
    src/main.cpp
)

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

add_executable(${TGT} ${FIGLET_SRC} ${SRC})
target_compile_features(${TGT} PRIVATE cxx_std_17)
if (MSVC)
  target_compile_features(${TGT} PRIVATE cxx_std_20)
  target_compile_definitions(${TGT} PRIVATE
    _USE_MATH_DEFINES = 1
  )
endif()
target_include_directories(${TGT} PRIVATE 
    ext/embedFiglet/src 
    ext/spdlog/include 
    ext/tclap-1.2.3/include 
    ext/libxdrfile/include
	  ext/inih
    ext/inih/cpp
    src
)
target_link_libraries(${TGT} PUBLIC 
    ${CMAKE_THREAD_LIBS_INIT}
    spdlog
    xdrfile
	  inih
)

if (USE_ONEAPI)
  target_compile_definitions(${TGT} PRIVATE USE_ONEAPI=1)
else()
  target_link_libraries(${TGT} PRIVATE OpenMP::OpenMP_CXX)
endif()

add_custom_command(
    TARGET ${TGT} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/${BIN}        
            ${CMAKE_SOURCE_DIR}/${BIN_OUT}
)

add_executable(jb_test ext/test_sims/jb_test.cpp)

if (Boost_FOUND)
  add_executable(ljmd_test 
      ext/test_sims/wblennardjones/src/coordinates.cpp 
      ext/test_sims/wblennardjones/src/main.cpp 
      ext/test_sims/wblennardjones/src/NeighborList.cpp 
      ext/test_sims/wblennardjones/src/PdbFile.cpp
      ext/test_sims/wblennardjones/src/Rdf.cpp 
      ext/test_sims/wblennardjones/src/ThermodynamicVariable.cpp 
      ext/test_sims/wblennardjones/src/Thermostat.cpp 
      ext/test_sims/wblennardjones/src/cubicbox.cpp 
      ext/test_sims/wblennardjones/src/utils.cpp
      ext/test_sims/wblennardjones/src/Velocity.cpp
  )
  target_include_directories(ljmd_test PRIVATE 
      ${BOOST_INCLUDES}
      ext/libxdrfile/include
      ext/test_sims/wblennardjones/include
  )
  target_link_libraries(ljmd_test PRIVATE 
    ${LIBBBOOST}
    xdrfile
  )
  if (!USE_ONEAPI)
    target_link_libraries(ljmd_test PRIVATE OpenMP::OpenMP_CXX)
  endif()
else()
  message(WARNING "Not building ljmd test sim since Boost was not found.")
endif()

MESSAGE( STATUS "C++ compiler                 = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE      = ${CMAKE_CXX_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_DEBUG        = ${CMAKE_CXX_FLAGS_DEBUG}" )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME            = ${CMAKE_SYSTEM_NAME}" )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR       = ${CMAKE_SYSTEM_PROCESSOR}" )
IF(APPLE)
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT          = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES    = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET= ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()
