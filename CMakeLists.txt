cmake_minimum_required(VERSION 3.14)
project(AlgorithmProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to choose CUDA or CPU. Default CPU
option(USE_CUDA "Use CUDA implementation" OFF)

# -----------------------------------------------------------------
# Smith-Waterman sources
set(ALG smith_waterman)

set(SRC
    algorithms/${ALG}/smith_waterman.cpp
    algorithms/${ALG}/cpu_sw_create_scoring_matrix.cpp
)

if(USE_CUDA)
    # Modern CMake way to find CUDA
    find_package(CUDAToolkit REQUIRED)
    list(APPEND SRC algorithms/${ALG}/cuda_sw_create_scoring_matrix.cu)
endif()

# Library for Smith-Waterman
add_library(${ALG}_lib ${SRC})

target_include_directories(${ALG}_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/algorithms/${ALG}
)

# Link CUDA runtime if needed
if(USE_CUDA)
    target_link_libraries(${ALG}_lib PRIVATE CUDA::cudart)
    target_compile_definitions(${ALG}_lib PUBLIC USE_CUDA=1)
endif()

# -----------------------------------------------------------------
# GoogleTest
find_package(GTest REQUIRED)

set(TEST_SRC
    tests/sw_test.cpp
)

add_executable(algorithm_tests ${TEST_SRC})

target_link_libraries(algorithm_tests PRIVATE ${ALG}_lib GTest::GTest GTest::Main)

enable_testing()
add_test(NAME AlgorithmTests COMMAND algorithm_tests)
