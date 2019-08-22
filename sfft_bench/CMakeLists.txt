cmake_minimum_required(VERSION 3.10)
project(sfft_bench)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

find_package(OpenMP REQUIRED)

set(SFFT_DIR "${CMAKE_SOURCE_DIR}/../")

include_directories(${SFFT_DIR}/src/)

add_executable(sfft_bench main.cpp)
target_link_libraries(sfft_bench OpenMP::OpenMP_CXX)
target_link_libraries(sfft_bench fftw3)
target_link_libraries(sfft_bench ${SFFT_DIR}/build/libsfft.so)

add_executable(measure measure.cpp)
target_link_libraries(measure OpenMP::OpenMP_CXX)
target_link_libraries(measure fftw3)
target_link_libraries(measure ${SFFT_DIR}/build/libsfft.so)
