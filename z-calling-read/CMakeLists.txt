cmake_minimum_required(VERSION 3.22)
project(z-calling-read)

set(CMAKE_CXX_STANDARD 17)

add_executable(z-calling-read main.cpp BamRead.h BamRead.cpp ../3rdparty/libsvm/svm.cpp func.cpp func.h)

link_directories(
        ../3rdparty/htslib/lib
)
include_directories(
        ../3rdparty/htslib/include
)


target_link_libraries(
        ${PROJECT_NAME}
        "libhts.so"
)
