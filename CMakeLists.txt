cmake_minimum_required(VERSION 2.8)

project(SlicedOptimalTransportSampler)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

ADD_LIBRARY(SOURCES
        Transport/NBallRadonManager.cpp
        Tools/mapping.cpp
        Tools/my_utility.cpp
        Math/myMath.cpp
        )

add_executable(sot Mains/slicedOptimNDBall.cpp)

find_package(GSL REQUIRED)
if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    target_link_libraries(SOURCES ${GSL_LIBRARIES})
endif()

find_package(OpenMP)
## On MACOS, you would need "libomp" (brew install libomp)
if(APPLE)
    MESSAGE(STATUS "Apple OpenMP build")
    SET(OpenMP_C_LIB_NAMES "omp")
    SET(OpenMP_CXX_LIB_NAMES "omp" )
    if(EXISTS "/usr/local/opt/libomp/lib/libomp.dylib")
        SET(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
        SET(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib")
        SET(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
    else()
        SET(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
        SET(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
        SET(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    endif()
endif()

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
SET(CMAKE_C_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")

target_link_libraries(sot SOURCES ${OpenMP_omp_LIBRARY})
