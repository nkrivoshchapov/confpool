cmake_minimum_required(VERSION 3.4...3.18)
project(confpool)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(pybind11)
add_subdirectory(fmt)

link_directories(fmt/build)
include_directories(fmt/include)

find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(GSL REQUIRED)
pybind11_add_module(confpool src/main.cpp)
target_link_libraries(confpool PRIVATE fmt::fmt GSL::gslcblas)
