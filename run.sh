#!/bin/bash
cmake ..
make
cd ./sfft_bench/cmake-build-debug/
make
./sfft_bench
