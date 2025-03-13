# CMake file for benchmark and validation of JetReconstruction.jl C-bindings
#
# (C) CERN, 2025
# MIT License

cmake_minimum_required(VERSION 3.12)
project(JetReconstructionBenchmarks)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(HepMC3 REQUIRED)
find_package(JetReconstruction REQUIRED)
find_package(Julia REQUIRED) # Find Julia for GC utilities

set(CMAKE_CXX_STANDARD 17)

# JetReconstruction.jl C-bindings executable runs reconstruction then outputs either
# exclusive or inclusive jets
add_executable(cjetreconstruction-finder
    src/cjetreconstruction-finder.cc
    src/cjetreconstruction-utils.cc
)

target_include_directories(cjetreconstruction-finder PRIVATE ${Julia_INCLUDE_DIRS})
target_link_libraries(cjetreconstruction-finder 
    HepMC3::HepMC3
    JetReconstruction::JetReconstruction ${Julia_LIBRARY}
)
