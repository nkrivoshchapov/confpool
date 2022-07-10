cmake_minimum_required(VERSION 3.4...3.18)
project(confpool)

add_subdirectory(pybind11)
pybind11_add_module(confpool src/main.cpp)
