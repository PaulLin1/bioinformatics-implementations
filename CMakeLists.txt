cmake_minimum_required(VERSION 3.14)
project(AlgorithmProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to choose CUDA or CPU. Default CPU
option(USE_CUDA "Use CUDA implementation" OFF)

# -----------------------------------------------------------------
# Hamming Distance library
set(HD_SRC
    src/algorithms/hamming_distance/hamming_distance.cpp
    src/algorithms/hamming_distance/cpu_hamming_distance.cpp
    )

if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    list(APPEND HD_SRC src/algorithms/hamming_distance/cuda_hamming_distance.cu)
endif()

add_library(hamming_distance_lib ${HD_SRC})
target_include_directories(hamming_distance_lib PUBLIC ${CMAKE_SOURCE_DIR}/include/)

if(USE_CUDA)
    target_link_libraries(hamming_distance_lib PRIVATE CUDA::cudart)
    target_compile_definitions(hamming_distance_lib PUBLIC USE_CUDA=1)
endif()

# -----------------------------------------------------------------
# Needleman Wunsch library
set(NW_SRC
    src/algorithms/needleman_wunsch/needleman_wunsch.cpp
    src/algorithms/needleman_wunsch/cpu_nw_create_scoring_matrix.cpp
    )

# if(USE_CUDA)
#     find_package(CUDAToolkit REQUIRED)
#     list(APPEND NW_SRC src/algorithms/needleman_wunsch/cuda_needleman_wunsch.cu)
# endif()

add_library(needleman_wunsch_lib ${NW_SRC})
target_include_directories(needleman_wunsch_lib PUBLIC ${CMAKE_SOURCE_DIR}/include/)

if(USE_CUDA)
    target_link_libraries(needleman_wunsch_lib PRIVATE CUDA::cudart)
    target_compile_definitions(needleman_wunsch_lib PUBLIC USE_CUDA=1)
endif()

# -----------------------------------------------------------------
# Smith-Waterman library
set(SW_SRC
    src/algorithms/smith_waterman/smith_waterman.cpp
    src/algorithms/smith_waterman/cpu_sw_create_scoring_matrix.cpp
)

if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    list(APPEND SW_SRC src/algorithms/smith_waterman/cuda_sw_create_scoring_matrix.cu)
endif()

add_library(smith_waterman_lib ${SW_SRC})
target_include_directories(smith_waterman_lib PUBLIC ${CMAKE_SOURCE_DIR}/include/)

if(USE_CUDA)
    target_link_libraries(smith_waterman_lib PRIVATE CUDA::cudart)
    target_compile_definitions(smith_waterman_lib PUBLIC USE_CUDA=1)
endif()

# -----------------------------------------------------------------
# GoogleTest for Hamming Distance
find_package(GTest REQUIRED)

set(HD_TEST_SRC tests/hd_test.cpp)

add_executable(hd_test ${HD_TEST_SRC})
target_link_libraries(hd_test PRIVATE hamming_distance_lib GTest::GTest GTest::Main)

enable_testing()
add_test(NAME HammingDistanceTests COMMAND hd_test)

# -----------------------------------------------------------------
# GoogleTest for needleman-Wunsch
set(NW_TEST_SRC tests/nw_test.cpp)

add_executable(nw_test ${NW_TEST_SRC})
target_link_libraries(nw_test PRIVATE needleman_wunsch_lib GTest::GTest GTest::Main)

enable_testing()
add_test(NAME needlemanWunschTests COMMAND nw_test)

# -----------------------------------------------------------------
# GoogleTest for Smith-Waterman
set(SW_TEST_SRC tests/sw_test.cpp)

add_executable(sw_test ${SW_TEST_SRC})
target_link_libraries(sw_test PRIVATE smith_waterman_lib GTest::GTest GTest::Main)

enable_testing()
add_test(NAME SmithWatermanTests COMMAND sw_test)
