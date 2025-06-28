cmake_minimum_required(VERSION 3.12)
project(BallisticsModel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

add_library(ballistics_lib STATIC
    src/coords.cpp
    src/utils.cpp
)

target_include_directories(ballistics_lib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


find_package(GSL REQUIRED)

add_executable(ballistics_tests
    src/test.cpp
)

add_executable(paris_solver
    src/paris_solver.cpp
)

target_link_libraries(ballistics_tests PRIVATE ballistics_lib)

target_link_libraries(paris_solver
    PRIVATE
    ballistics_lib
    GSL::gsl
    GSL::gslcblas
)