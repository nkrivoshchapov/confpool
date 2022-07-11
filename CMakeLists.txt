cmake_minimum_required(VERSION 3.4...3.18)
project(confpool)

add_subdirectory(pybind11)
FIND_PACKAGE(Boost COMPONENTS system filesystem REQUIRED)

pybind11_add_module(confpool src/main.cpp)
