cmake_minimum_required(VERSION 3.16)
project(filtration_solver LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-DTRILIBRARY)

include_directories(external)
include_directories(external/Eigen)

add_executable(solver
    src/main.cpp
    external/triangle/triangle.c
)

target_include_directories(solver PRIVATE
    external/triangle
    external
)

if(WIN32)
    target_compile_definitions(solver PRIVATE NO_TIMER)
endif()