cmake_minimum_required(VERSION 3.10)
project(lidDrivenCavity)

# C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include headers
include_directories(include)

# Gather source files
file(GLOB SRC_FILES src/*.cpp lidDrivenCavity.cpp)
# Executable
add_executable(lidDrivenCavity ${SRC_FILES})

