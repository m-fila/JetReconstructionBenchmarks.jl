# CMake file for fastjet codes used for benchmarks and validation
#
# (C) CERN, 2024-2025
# MIT License

cmake_minimum_required(VERSION 3.12)
project(FastJetBenchmarks)

find_package(HepMC3 REQUIRED)
find_package(fastjet REQUIRED)
find_package(siscone REQUIRED)

set(CMAKE_CXX_STANDARD 17)

# Fastjet executable runs reconstruction then outputs either
# exclusive or inclusive jets
add_executable(fastjet-finder
    src/fastjet-finder.cc
    src/fastjet-utils.cc
)

target_link_libraries(fastjet-finder 
    HepMC3::HepMC3 fastjet::fastjet
)
