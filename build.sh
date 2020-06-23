#!/bin/sh

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_TESTING=ON -DCMAKE_CXX_COMPILER=clang++ ..
make