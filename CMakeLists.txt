cmake_minimum_required(VERSION 3.14)
project(AlgorithmProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to choose CUDA or CPU. Default CPU
option(USE_CUDA "Use CUDA implementation" OFF)
if(USE_CUDA)
    enable_language(CUDA)
endif()

# List of algorithms
set(ALGORITHMS
    smith_waterman
    # needleman_wunsch
    # edit_distance
)

foreach(ALG ${ALGORITHMS})
    if(USE_CUDA)
        set(SRC algorithm/${ALG}/cuda_${ALG}.cu)
    else()
        set(SRC algorithm/${ALG}/naive_${ALG}.cpp)
    endif()

    # Main executable for each algorithm
    add_executable(${ALG}_program main.cpp ${SRC})

    # Library for each algorithm (optional)
    add_library(${ALG}_lib
        algorithm/${ALG}/${ALG}.cpp
        utils/utils.cpp
    )
    target_include_directories(${ALG}_lib PUBLIC
        ${CMAKE_SOURCE_DIR}/algorithm/${ALG}
        ${CMAKE_SOURCE_DIR}/utils
    )

    target_link_libraries(${ALG}_program PRIVATE ${ALG}_lib)
endforeach()

# GoogleTest
# find_package(GTest REQUIRED)
# add_executable(algorithm_tests tests/algorithm_tests.cpp)
# foreach(ALG ${ALGORITHMS})
#    target_link_libraries(algorithm_tests PRIVATE ${ALG}_lib)
# endforeach()

enable_testing()
add_test(NAME AlgorithmTests COMMAND algorithm_tests)
