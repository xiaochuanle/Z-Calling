cmake_minimum_required(VERSION 3.22)
project(z-seq)

set(CMAKE_CXX_STANDARD 17)

add_executable(z-seq main.cpp)

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
