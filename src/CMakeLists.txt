# Set include directory.
include_directories("${CMAKE_SOURCE_DIR}/lib")

# Add compiler switch to generate SSE2 instructions.
# Note: This should be moved to CMakeLists.txt in the root directory.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")

# Create executable: mmserve.
add_executable(mmserve mmserve.cpp)
target_link_libraries(mmserve mastermind)

# Create executable: mmstrat.
add_executable(mmstrat mmstrat.cpp)
target_link_libraries(mmstrat mastermind)

