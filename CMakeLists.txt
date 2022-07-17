cmake_minimum_required(VERSION 3.4...3.18)
project(confpool)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory(pybind11)
add_subdirectory(fmt)
add_subdirectory(eigen)

link_directories(fmt/build)
include_directories(fmt/include)
include_directories(eigen/Eigen)

find_package (Eigen3 REQUIRED NO_MODULE)
find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(GSL REQUIRED)
pybind11_add_module(confpool src/main.cpp)
target_link_libraries(confpool PRIVATE fmt::fmt Eigen3::Eigen GSL::gslcblas)
