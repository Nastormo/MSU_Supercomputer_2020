#!/bin/bash
echo "Nvcc"
nvcc -ccbin mpicxx -std=c++11 block.cu process.cpp function.cpp config.cpp cuda.cu main.cpp -o cuda_main
echo "Done"
echo "Omp"
mpicxx -openmp block.cpp process.cpp function.cpp config.cpp omp.cpp main.cpp -o omp_main
echo "Done"
echo "None"
mpicxx block.cpp process.cpp function.cpp config.cpp omp.cpp main.cpp -o main
echo "Done"