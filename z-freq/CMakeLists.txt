cmake_minimum_required(VERSION 3.22)
project(zfreq)

set(CMAKE_CXX_STANDARD 17)

add_executable(zfreq zfreq.cpp)

target_link_libraries(
        ${PROJECT_NAME}
)
