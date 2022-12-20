add_custom_target(test 
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/test-mmstrat.pl"
    "${CMAKE_BINARY_DIR}/bin/mmstrat"
  DEPENDS mmstrat)
