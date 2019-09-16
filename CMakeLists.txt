cmake_minimum_required(VERSION 3.1)

project(mincv)


set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

find_path(MPFR_INCLUDE_DIR mpfr.h PATHS /include /usr/include /usr/local/include)
find_library(MPFR_LIBRARY NAMES mpfr PATHS /lib /usr/lib /usr/local/lib)

include_directories(${CMAKE_SOURCE_DIR})
include_directories("/usr/include/eigen3")
include_directories("/usr/local/include/libcmaes")
include_directories(${MPFR_INCLUDE_DIR})
link_directories("/usr/local/lib")

add_executable(${PROJECT_NAME} "main.cpp")
target_link_libraries(${PROJECT_NAME} ${MPFR_LIBRARY} cmaes)
