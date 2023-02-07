cmake_minimum_required(VERSION 3.24)
project(Mastermind)

# Set C++ version.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Treat warnings as errors.
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

# Find source files.
file(GLOB SOURCES "*.cpp")

# Create executable.
add_executable(mastermind ${SOURCES})

# Turn on warnings.
target_compile_options(mastermind PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
)
