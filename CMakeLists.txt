cmake_minimum_required(VERSION 3.22)

## project name
project(Z-Calling)

set(CMAKE_CXX_STANDARD 17)

MACRO(SUBDIRLIST result curdir)
    FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
    SET(dirlist "")
    FOREACH (child ${children})
        IF (IS_DIRECTORY ${curdir}/${child})
            LIST(APPEND dirlist ${child})
        ENDIF ()
    ENDFOREACH ()
    SET(${result} ${dirlist})
ENDMACRO()

find_package(Torch REQUIRED)
find_package(ZLIB REQUIRED)

link_directories(
        3rdparty/htslib/lib
)
include_directories(
        3rdparty/htslib/include
)


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
        3rdparty/libsvm/svm.cpp
        z-calling-read/func.cpp
        z-calling-read/func.h)


SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

foreach (dirname IN LISTS SUBDIRS)
    if (${dirname} MATCHES "z-bam2txt")
        add_executable(z-bam2txt ${z-bam2txt_LIST})
        target_link_libraries(z-bam2txt "libhts.so" ${ZLIB_LIBRARIES})
    elseif (${dirname} MATCHES "z-calling-base")
        add_executable(z-calling-base ${z-calling-base_LIST})
        target_link_libraries(z-calling-base "libhts.so" ${TORCH_LIBRARIES} ${ZLIB_LIBRARIES})
    elseif (${dirname} MATCHES "z-calling-read")
        add_executable(z-calling-read ${z-calling-read_LIST})
        target_link_libraries(z-calling-read "libhts.so")
    elseif (${dirname} MATCHES "z-freq")
        add_executable(zfreq z-freq/zfreq.cpp)
        target_link_libraries(zfreq)
    elseif (${dirname} MATCHES "z-seq")
        add_executable(z-seq z-seq/main.cpp)
        target_link_libraries(z-seq "libhts.so")
    endif ()
endforeach ()
