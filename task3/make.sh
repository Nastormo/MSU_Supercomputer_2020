#!/bin/bash
#mpic++ -std=c++98 block.cpp process.cpp function.cpp config.cpp main.cpp -o main
nvcc -ccbin mpicxx -std=c++11 block.cpp process.cpp function.cpp config.cpp main.cpp -o cuda_main
echo "Done"