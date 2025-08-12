cmake_minimum_required(VERSION 3.14)
project(SmithWatermanProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add GoogleTest
find_package(GTest REQUIRED)

# Library target
add_library(smith_waterman smith_waterman/smith_waterman.cpp)
target_include_directories(smith_waterman PUBLIC ${CMAKE_SOURCE_DIR}/smith_waterman)

# Test target
add_executable(smith_waterman_test tests/smith_waterman_test.cpp)
target_link_libraries(smith_waterman_test PRIVATE smith_waterman GTest::gtest_main)

# Enable tests
enable_testing()
add_test(NAME SmithWatermanTests COMMAND smith_waterman_test)
