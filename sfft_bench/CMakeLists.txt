cmake_minimum_required(VERSION 3.12)
project(sfft_bench)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)

add_executable(sfft_bench main.cpp)
target_link_libraries(sfft_bench OpenMP::OpenMP_CXX)
target_link_libraries(sfft_bench fftw3)

set(SFFT_DIR "${CMAKE_SOURCE_DIR}/../")

include_directories(${SFFT_DIR}/src/)
target_link_libraries(sfft_bench ${SFFT_DIR}/build/libsfft.so)