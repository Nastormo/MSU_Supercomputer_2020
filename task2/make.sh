#!/bin/bash
mpic++ -std=c++98 -openmp block.cpp process.cpp function.cpp main.cpp -o main
#mpic++ -std=c++98 block.cpp process.cpp function.cpp main.cpp -o main
echo "Done"