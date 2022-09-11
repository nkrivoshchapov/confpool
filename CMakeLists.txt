cmake_minimum_required(VERSION 3.4...3.18)
project(confpool)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set (Boost_USE_STATIC_LIBS OFF CACHE BOOL "use static libraries from Boost")

add_subdirectory(pybind11)
add_subdirectory(fmt)
add_subdirectory(eigen)

link_directories(fmt/build)
include_directories(fmt/include)
include_directories(eigen/Eigen)
# add_compile_options("-W" "-Wall" "-O3" "-pedantic")


if (WIN32)
# disable autolinking in boost
add_definitions( -DBOOST_ALL_NO_LIB )

# force all boost libraries to dynamic link (we already disabled
# autolinking, so I don't know why we need this, but we do!)
add_definitions( -DBOOST_ALL_DYN_LINK )
endif()

find_package (Eigen3 REQUIRED NO_MODULE)
find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(GSL REQUIRED)
pybind11_add_module(confpool src/main.cpp src/confpool.cpp src/utils.cpp src/molproxy.cpp src/hungarian.cpp)
target_compile_options(confpool PRIVATE -pedantic -O3 -D NDEBUG -DBOOST_UBLAS_NDEBUG)
target_link_options(confpool PRIVATE -pedantic -O3 -D NDEBUG -DBOOST_UBLAS_NDEBUG)
target_link_libraries(confpool PRIVATE fmt::fmt Eigen3::Eigen GSL::gslcblas)
