cmake_minimum_required(VERSION 3.22)

project(Z-Calling)

set(CMAKE_CXX_STANDARD 17)

# ==========================================
# 1. AUTOMATIC DEPENDENCY DETECTION (NEW)
# ==========================================

# A. Find PyTorch (Dynamically from Python)
execute_process(
    COMMAND python3 -c "import torch; print(torch.utils.cmake_prefix_path)"
    OUTPUT_VARIABLE TORCH_CMAKE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
list(APPEND CMAKE_PREFIX_PATH "${TORCH_CMAKE_PATH}")
find_package(Torch REQUIRED)

# B. Find HTSlib (From Conda environment)
find_path(HTSLIB_INCLUDE_DIR htslib/hts.h)
find_library(HTSLIB_LIBRARY hts)
if(NOT HTSLIB_LIBRARY)
    message(FATAL_ERROR "htslib not found! Did you activate the conda environment?")
endif()

# C. Standard Libraries
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Global Include Directories
include_directories(${HTSLIB_INCLUDE_DIR})

# ==========================================
# 2. SOURCE FILE DEFINITIONS (KEPT FROM YOUR ORIGINAL)
# ==========================================

SET(COMMON_LIST
        corelib/fasta.cpp
        corelib/hbn_aux.c
        corelib/line_reader.cpp
        corelib/sam_parser.cpp
        corelib/split_string_by_char.cpp
        corelib/getMemorySize.c
        corelib/getRSS.c
        str_util/c_ncbi_blast_aux.c
        str_util/ncbi_blast_aux.cpp
        str_util/ncbistr_util.cpp
        str_util/numeric_str_interconv.cpp
        str_util/str_cmp.cpp
        str_util/str_util.cpp)

SET(z-bam2txt_LIST
        ${COMMON_LIST}
        z-bam2txt/bam_to_txt.cpp
        z-bam2txt/necat_info.cpp
        z-bam2txt/sam_map_info.cpp)

SET(z-calling-base_LIST
        ${COMMON_LIST}
        z-calling-base/main.cpp
        z-calling-base/build_mod_bam.cpp
        z-calling-base/make_read_features.cpp
        z-calling-base/necat_info.cpp)

SET(z-calling-read_LIST
        z-calling-read/main.cpp
        z-calling-read/BamRead.h
        z-calling-read/BamRead.cpp
        3rdparty/libsvm/svm.cpp  # Keeping your bundled SVM to avoid code changes
        z-calling-read/func.cpp
        z-calling-read/func.h)

# ==========================================
# 3. TARGET DEFINITIONS (UPDATED LINKING)
# ==========================================

# Target: z-bam2txt
add_executable(z-bam2txt ${z-bam2txt_LIST})
target_link_libraries(z-bam2txt ${HTSLIB_LIBRARY} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Target: z-calling-base (Needs Torch)
add_executable(z-calling-base ${z-calling-base_LIST})
target_link_libraries(z-calling-base ${HTSLIB_LIBRARY} "${TORCH_LIBRARIES}" ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Target: z-calling-read
add_executable(z-calling-read ${z-calling-read_LIST})
target_link_libraries(z-calling-read ${HTSLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

# Target: zfreq
# Check if directory exists to avoid errors
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/z-freq")
    add_executable(zfreq z-freq/zfreq.cpp)
    target_link_libraries(zfreq ${HTSLIB_LIBRARY}) # Assuming it needs HTSlib, removed if not
endif()

# Target: z-seq
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/z-seq")
    add_executable(z-seq z-seq/main.cpp)
    target_link_libraries(z-seq ${HTSLIB_LIBRARY})
endif()
