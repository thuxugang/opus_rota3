cmake_minimum_required(VERSION 3.5)
project(rota3)

set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O2")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -O2")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

message("=================================")
message(${EIGEN3_INCLUDE_DIR})

aux_source_directory(. DIR_SRCS)
add_executable(rota3 ${DIR_SRCS})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(hiredis REQUIRED)
include_directories( ${HIREDIS_INCLUDE_DIRS} )
target_link_libraries( rota3 ${HIREDIS_LIBRARIES} )

message("=================================")
message(${HIREDIS_INCLUDE_DIRS})