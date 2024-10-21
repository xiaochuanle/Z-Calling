cmake_minimum_required(VERSION 3.22)
project(z-bam2txt)

set(CMAKE_CXX_STANDARD 17)

find_package(ZLIB REQUIRED)

add_executable(z-bam2txt bam_to_txt.cpp necat_info.cpp sam_map_info.cpp
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
        ${ZLIB_LIBRARIES}
)
