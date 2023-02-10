cmake_minimum_required(VERSION 3.24)
project(Mastermind)

# Set C++ version.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

# Add tests.
include(CTest)
add_test(NAME test_count_6x4  COMMAND $<TARGET_FILE:mastermind> -r 6x4 count)
add_test(NAME test_count_10p4 COMMAND $<TARGET_FILE:mastermind> -r 10p4 count)
add_test(NAME test_list_6x4   COMMAND $<TARGET_FILE:mastermind> -r 6x4 list)
add_test(NAME test_list_10p4  COMMAND $<TARGET_FILE:mastermind> -r 10p4 list)
add_test(NAME test_filter     COMMAND $<TARGET_FILE:mastermind> -c 0123:2A0B -c 1234:1A1B list)

