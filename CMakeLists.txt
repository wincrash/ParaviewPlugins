project(ParaViewPlugins) 
cmake_minimum_required(VERSION 3.9)
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory("HDF5Reader")
add_subdirectory("Voronoi")
add_subdirectory("CrackAlgorithms")
