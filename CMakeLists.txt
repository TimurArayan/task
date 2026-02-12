cmake_minimum_required(VERSION 3.16)
project(filtration_solver LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(solver
    src/main.cpp
    src/SymmetricSparseMatrices.cpp
    external/triangle/triangle.c
)

target_include_directories(solver PRIVATE
    external/triangle
    external
)

target_compile_definitions(solver PRIVATE
    TRILIBRARY
    ANSI_DECLARATORS
    NO_TIMER
)
