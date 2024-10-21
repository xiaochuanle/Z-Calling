cmake_minimum_required(VERSION 3.22)
project(z-calling-base)

set(CMAKE_CXX_STANDARD 17)

set(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda-11)
set(CUDNN_INCLUDE_DIR /usr/local/cuda-11/include)
set(CUDNN_LIBRARY /usr/local/cuda-11/lib64)
set(CMAKE_CUDA_COMPILER /usr/local/cuda-11.8/bin/nvcc)
set(CMAKE_PREFIX_PATH /home/username/Software/libtorch/share/cmake)
set(CMAKE_PREFIX_PATH /usr/lib/x86_64-linux-gnu)
set(CAFFE2_USE_CUDNN TRUE)
set(Torch_DIR /home/username/Software/libtorch/share/cmake/Torch)
find_package(Torch REQUIRED)
find_package(ZLIB REQUIRED)

add_executable(z-calling-base main.cpp build_mod_bam.cpp make_read_features.cpp necat_info.cpp
        ../corelib/fasta.cpp
        ../corelib/hbn_aux.c
        ../corelib/line_reader.cpp
        ../corelib/sam_parser.cpp
        ../corelib/split_string_by_char.cpp
        ../corelib/getMemorySize.c
        ../corelib/getRSS.c
        ../str_util/c_ncbi_blast_aux.c
        ../str_util/ncbi_blast_aux.cpp
        ../str_util/ncbistr_util.cpp
        ../str_util/numeric_str_interconv.cpp
        ../str_util/str_cmp.cpp
        ../str_util/str_util.cpp)

link_directories(
        ../3rdparty/htslib/lib
)
include_directories(
        ../3rdparty/htslib/include
)

target_link_libraries(
        ${PROJECT_NAME}
        "libhts.so"
        ${TORCH_LIBRARIES}
        ${ZLIB_LIBRARIES}
)
