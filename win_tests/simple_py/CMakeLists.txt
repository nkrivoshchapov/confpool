cmake_minimum_required(VERSION 3.4...3.18)
project(confpool)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
find_package (pybind11 REQUIRED)
pybind11_add_module(confpool main.cpp)
