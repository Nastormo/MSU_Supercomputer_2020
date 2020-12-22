#!/bin/bash
#mpic++ -std=c++98 block.cpp process.cpp function.cpp config.cpp main.cpp -o main
echo "Nvcc"
nvcc -ccbin mpicxx -std=c++11 block.cpp process.cpp function.cpp config.cpp cuda.cu main.cpp -o cuda_main
echo "Done"
echo "Omp"
mpicxx -std=c++11 -openmp block.cpp process.cpp function.cpp config.cpp omp.cpp main.cpp -o omp_main
echo "Done"