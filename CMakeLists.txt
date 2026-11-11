cmake_minimum_required(VERSION 3.15)
project(EMTSim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Solver library (everything except the circuit definition and main)
add_library(emtsim_core STATIC
  src/linalg.cpp
  src/elements.cpp
  src/network.cpp
  src/netlist.cpp
  src/switching.cpp
)
target_include_directories(emtsim_core PUBLIC src)
if(MSVC)
  target_compile_options(emtsim_core PUBLIC /W4 /permissive- /EHsc)
else()
  target_compile_options(emtsim_core PUBLIC -Wall -Wextra -Wpedantic)
endif()

# Simulator: the circuit in src/circuit.cpp
add_executable(emtsim src/main.cpp src/circuit.cpp)
target_link_libraries(emtsim PRIVATE emtsim_core)

# Tests against analytical solutions
enable_testing()
add_executable(test_emtsim tests/test_emtsim.cpp)
target_link_libraries(test_emtsim PRIVATE emtsim_core)
add_test(NAME analytical COMMAND test_emtsim)
