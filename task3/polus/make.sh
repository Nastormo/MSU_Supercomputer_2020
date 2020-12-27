module load SpectrumMPI
module load OpenMPI
echo "Nvcc"
nvcc -ccbin mpixlC -std=c++11 block.cu process.cpp function.cpp config.cpp cuda.cu main.cpp -o cuda_main
echo "Done"
echo "Omp"
mpixlC -std=c++11 -qsmp=omp block.cpp process.cpp function.cpp config.cpp omp.cpp main.cpp -o omp_main
echo "Done"
echo "None"
mpixlC -std=c++11 block.cpp process.cpp function.cpp config.cpp omp.cpp main.cpp -o main
echo "Done"