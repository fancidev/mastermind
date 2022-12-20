# Check for SSE2 support on the platform.
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(SSE 2.0 REQUIRED)
#if(NOT SSE_FOUND)
#  message(FATAL_ERROR "SSE support is missing on this architecture.")
#else()
#  message(STATUS "Found SSE intrinsic support: version ${SSE_VERSION}")
#endif()

# Add compiler switch to generate SSE2 instructions.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")

# List of source files.
set(SRC_LIST CodeBreaker.cpp Engine.cpp ObviousStrategy.cpp Codeword.cpp OptimalCodeBreaker.cpp ColorEquivalence.cpp Generation.cpp StrategyTree.cpp Compare.cpp ConstraintEquivalence.cpp DummyEquivalenceFilter.cpp Mask.cpp)

# Create static library.
add_library(mastermind STATIC ${SRC_LIST})
