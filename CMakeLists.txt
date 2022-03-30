project(SlicedOptimalTransportSampler)

cmake_minimum_required(VERSION 3.12)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

include(openmp)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})


ADD_LIBRARY(SOURCES
        Transport/NBallRadonManager.cpp
        Tools/mapping.cpp
        Tools/my_utility.cpp
        Math/myMath.cpp
        )
target_link_libraries(SOURCES GSL::gsl)

add_executable(sot Mains/slicedOptimNDBall.cpp)
target_link_libraries(sot SOURCES OpenMP::OpenMP_CXX)
