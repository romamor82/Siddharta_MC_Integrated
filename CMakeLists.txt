################################################################################
## CMake script for Siddharta Geant4 Monte Carlo simulation software
################################################################################

cmake_minimum_required (VERSION 2.6)
project(sidd_mc)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

################################################################################
## Specify source folders
set(FOLDERS_WITH_SOURCE
  src
  include
)
################################################################################
## Provide additional files to copy
set(SCRIPT_FILES
  CARD.DAT
  geom.txt
  touschek.dat
  vis.mac
  cfg_D_1_4dens_475deg.dat
  run.sh
)

################################################################################
## Point SOURCES and headers
foreach(dir ${FOLDERS_WITH_SOURCE})
  file(GLOB HEADERS_TMP ${dir}/*.h)
  file(GLOB SOURCES_TMP ${dir}/*.cpp)
  list(APPEND HEADERS ${HEADERS_TMP})
  list(APPEND SOURCES ${SOURCES_TMP})
endforeach()

list(APPEND FOLDERS ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_LIBRARY_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_CXX_FLAGS "-std=c++11 -Wunused-parameter ")

################################################################################
## Project version number
set(VERSION_MAJOR 3 CACHE STRING "Project major version number.")
set(VERSION_MINOR 0 CACHE STRING "Project minor version number.")
set(VERSION_PATCH 0 CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

################################################################################
## Include GEANT4 WITH VIS DRIVERS
## remember to compile geant4 with following flags
## -DGEANT4_USE_QT=ON -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_GDML=ON
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

################################################################################
## Include ROOT
find_package(ROOT REQUIRED)
message(STATUS "ROOT version: ${ROOT_VERSION}")

set(ROOT_CXX_FLAGS "-std=c++11 -Wunused-parameter")
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(SYSTEM ${Geant4_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})

################################################################################
## Create a main program using the library
add_executable(sidd_mc Siddharta.cpp ${SOURCES})
target_link_libraries(
  sidd_mc
  ${Geant4_LIBRARIES}
  ${ROOT_LIBRARIES}
)

## Copy additional files to bin directory
foreach(file_i ${SCRIPT_FILES})
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${file_i}
    ${CMAKE_CURRENT_BINARY_DIR}/bin/${file_i}
    COPYONLY
  )
endforeach(file_i)
