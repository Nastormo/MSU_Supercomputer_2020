#include "parallel.hpp"

#include <vector>
#include <iostream>
#include <math.h>

#include <cuda.h>

#define checkCudaErrors(val) check( (val), #val, __FILE__, __LINE__)

template<typename T>
void check(T err, const char* const func, const char* const file, const int line) {
  if (err != cudaSuccess) {
    std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
    std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
    exit(1);
  }
}

__global__ 
void cuda_test(int* vec)
{
  int idx = threadIdx.x;
  vec[idx] = vec[idx] + 1;
}

void test() {
    int size = 128;
    std::vector<int> vec (size, 0);

    int* devVec;

    checkCudaErrors(cudaMalloc((void**)&devVec, sizeof(int) * size));

    cudaMemcpy(devVec, vec.data(), sizeof(int) * size, cudaMemcpyHostToDevice);

    dim3 gridSize = dim3(1, 1, 1);
    dim3 blockSize = dim3(size, 1, 1);
    cuda_test<<<gridSize, blockSize>>>(devVec);

    cudaEvent_t syncEvent;

    cudaEventCreate(&syncEvent);
    cudaEventRecord(syncEvent, 0);
    cudaEventSynchronize(syncEvent);
    cudaMemcpy(vec.data(), devVec, sizeof(int) * size, cudaMemcpyDeviceToHost);

    std::cout << vec[10] << std::endl;
}

__device__
double u(double Lx, double Ly, double Lz, double x, double y, double z) {
    return sin((M_PI / Lx) * x) * 
        sin((M_PI / Ly) * y) * 
        sin((M_PI / Lz) * z);
}

__global__
void u0(double* block, double Lx, double Ly, double Lz, 
    int minI, int minJ, int minK, 
    double shiftX, double shiftY, double shiftZ) 
{
    
}

void init_u0(Block &b, Function3D &u) {
    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();

    int sizeBI = 1;
    int sizeBJ = 1;
    int sizeBK = sizeK - 2;

    checkCudaErrors(cudaMalloc((void**)&devVec, sizeof(double) * size));

    cudaMemcpy(devVec, vec.data(), sizeof(int) * size, cudaMemcpyHostToDevice);

    dim3 grid = dim3(sizeI / sizeBI, sizeJ / sizeBJ, sizeK / sizeBK);
    dim3 block = dim3(sizeBI, sizeBJ, sizeBK);
    
    std::cout << grid.x << grid.y << grid.z << std::endl;
    std::cout << block.x << block.y << block.z << std::endl;

    cuda_test<<<gridSize, blockSize>>>(devVec);

    cudaEvent_t syncEvent;

    cudaEventCreate(&syncEvent);
    cudaEventRecord(syncEvent, 0);
    cudaEventSynchronize(syncEvent);
    cudaMemcpy(vec.data(), devVec, sizeof(int) * size, cudaMemcpyDeviceToHost);



    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();

    #pragma omp parallel for
    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
                b.getElem(i, j, k) = u(b.getX(i), b.getY(j), b.getZ(k), 0);
            }
        }
    }
}

void init_u1(Block &b, const Block &u0, double tau, Function3D &u) {
    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();

    #pragma omp parallel for
    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
                b.getElem(i, j, k) = u0.getValElem(i, j, k) + (pow(tau, 2) / 2) * u0.lap_h(i, j, k);
            }
        }
    }
}

void step(Block &u2, const Block& u1, const Block& u0, double tau, Function3D &u) {
    int sizeI = u2.getSizeI();
    int sizeJ = u2.getSizeJ();
    int sizeK = u2.getSizeK();

    #pragma omp parallel for
    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
                u2.getElem(i, j, k) = 2 * u1.getValElem(i, j, k) - u0.getValElem(i, j, k) + pow(tau, 2) * u1.lap_h(i, j, k);
            }
        }
    }
}