cmake_minimum_required(VERSION 3.10)
project(Mastermind)

# Set C++ version.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler flags.
#target_compile_options(mastermind INTERFACE "-Wall")
# Add compiler switch.
# Note: This should be moved to CMakeLists.txt in the root directory.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Create executable: mmserve.
add_executable(mastermind canonical.cpp codemaker.cpp main.cpp
               codebreaker.cpp heuristicbreaker.cpp)
