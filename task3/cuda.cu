#include "parallel.hpp"

#include <vector>
#include <iostream>
#include <math.h>

#include <cuda.h>

#define SAFE_CALL( CallInstruction ) { \
    cudaError_t cuerr = CallInstruction; \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error: %s at call \"" #CallInstruction "\"\n", cudaGetErrorString(cuerr)); \
            throw "error in CUDA API function, aborting..."; \
    } \
}

#define SAFE_KERNEL_CALL( KernelCallInstruction ) { \
    KernelCallInstruction; \
    cudaError_t cuerr = cudaGetLastError(); \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error in kernel launch: %s at kernel \"" #KernelCallInstruction "\"\n", cudaGetErrorString(cuerr)); \
            throw "error in CUDA kernel launch, aborting..."; \
    } \
    cuerr = cudaDeviceSynchronize(); \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error in kernel execution: %s at kernel \"" #KernelCallInstruction "\"\n", cudaGetErrorString(cuerr)); \
            throw "error in CUDA kernel execution, aborting..."; \
    } \
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

    SAFE_CALL(cudaMalloc((void**)&devVec, sizeof(int) * size));

    SAFE_CALL(cudaMemcpy(devVec, vec.data(), sizeof(int) * size, cudaMemcpyHostToDevice));

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
double u(double Lx, double Ly, double Lz, double a, double x, double y, double z, double t) {
    return sin((M_PI / Lx) * x) * 
        sin((M_PI / Ly) * y) * 
        sin((M_PI / Lz) * z) *
        cos(a * t);;
}

__global__
void u0(double* block, 
    double Lx, double Ly, double Lz, double a,
    double shiftX, double shiftY, double shiftZ,
    int sizeI, int sizeJ, int sizeK,
    int minI, int minJ, int minK)
{
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    block[i * (sizeJ * sizeK) + j * (sizeK) + k] = u(Lx, Lz, Lz, a,
        (i + minI) * shiftX, (j + minJ) * shiftY, (k + minK) * shiftZ, 0.0f);
}

void init_u0(Block &b, Function3D &u) {
    double Lx = u.getLx();
    double Ly = u.getLy();
    double Lz = u.getLz();
    double a = u.a_t();

    double shiftX = b.getShiftX();
    double shiftY = b.getShiftY();
    double shiftZ = b.getShiftZ();

    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();

    int minI = b.getMinI();
    int minJ = b.getMinJ();
    int minK = b.getMinK();

    int size = sizeI * sizeJ * sizeK;

    int sizeBI = 1;
    int sizeBJ = 1;
    int sizeBK = sizeK - 2;

    double* d_block;

    SAFE_CALL(cudaMalloc((void**)&d_block, sizeof(double) * size));

    dim3 grid = dim3((sizeK - 2) / sizeBK, (sizeJ - 2) / sizeBJ, (sizeI - 2) / sizeBI);
    dim3 block = dim3(sizeBK, sizeBJ, sizeBI);

    std::cout << grid.x << ' ' << grid.y << ' ' << grid.z << std::endl;
    std::cout << block.x << ' ' << block.y << ' ' << block.z << std::endl;

    u0<<<grid, block>>>(d_block, 
        Lx, Ly, Lz, a,
        shiftX, shiftY, shiftZ, 
        sizeI, sizeJ, sizeK, 
        minI, minJ, minK);

    SAFE_CALL(cudaMemcpy(b.getData().data(), d_block, sizeof(double) * size, cudaMemcpyDeviceToHost));
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