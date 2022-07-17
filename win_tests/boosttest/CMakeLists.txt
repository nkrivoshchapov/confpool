cmake_minimum_required(VERSION 3.4...3.18)
project(confpool)

find_package(GSL REQUIRED)
add_executable(confpool main.cpp)
target_link_libraries(confpool PRIVATE GSL::gslcblas)
