cmake_minimum_required(VERSION 3.14)

project(AlgorithmProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_CUDA "Use CUDA implementation" OFF)

if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    add_compile_definitions(USE_CUDA=1)
    set(CMAKE_CUDA_ARCHITECTURES 90)
    project(AlgorithmProject LANGUAGES CXX CUDA)

endif()

if(USE_CUDA)
    set(CMAKE_CUDA_ARCHITECTURES 90)
endif()


# ----------------- Hamming Distance -----------------
set(HD_SRC
    src/algorithms/hamming_distance/cpu_hamming_distance.cpp
)

add_library(hamming_distance_lib ${HD_SRC})
target_include_directories(hamming_distance_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/algorithms/hamming_distance
)

# ----------------- Needleman-Wunsch -----------------
set(NW_SRC
    src/algorithms/needleman_wunsch/cpu_needleman_wunsch.cpp
)

add_library(needleman_wunsch_lib ${NW_SRC})
target_include_directories(needleman_wunsch_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/algorithms/needleman_wunsch
)

# ----------------- Smith-Waterman -----------------
set(SW_SRC
    src/algorithms/smith_waterman/cpu_smith_waterman.cpp
)

add_library(smith_waterman_lib ${SW_SRC})
target_include_directories(smith_waterman_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/algorithms/smith_waterman
)

# ----------------- Executable -----------------
set(EXE_SRC src/pipeline/pipeline.cpp)

if(USE_CUDA)
    list(APPEND EXE_SRC
        src/algorithms/hamming_distance/cuda_hamming_distance.cu
        src/algorithms/needleman_wunsch/cuda_needleman_wunsch.cu
        src/algorithms/smith_waterman/cuda_smith_waterman.cu
    )
endif()

add_executable(pipeline ${EXE_SRC})
target_include_directories(pipeline PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(pipeline
    PRIVATE
        hamming_distance_lib
        needleman_wunsch_lib
        smith_waterman_lib
        $<$<BOOL:${USE_CUDA}>:CUDA::cudart>
)
