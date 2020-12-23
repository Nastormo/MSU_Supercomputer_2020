#include "parallel.hpp"

#include <vector>
#include <iostream>
#include <math.h>

#include <cuda.h>

#include <mpi.h>

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

__constant__ double Lx;
__constant__ double Ly;
__constant__ double Lz;
__constant__ double a_t;

__constant__ double shiftX;
__constant__ double shiftY;
__constant__ double shiftZ;

__constant__ int sizeI;
__constant__ int sizeJ;
__constant__ int sizeK;

__constant__ int minI;
__constant__ int minJ;
__constant__ int minK;

__device__
double u(double x, double y, double z, double t) {
    return sin((M_PI / Lx) * x) * 
           sin((M_PI / Ly) * y) * 
           sin((M_PI / Lz) * z) *
           cos(a_t * t);;
}

__global__
void u0(double* block)
{
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    block[i * (sizeJ * sizeK) + j * (sizeK) + k] = u((i + minI) * shiftX, (j + minJ) * shiftY, (k + minK) * shiftZ, 0.0f);
}

void init_u0(Block &b, Function3D &u) {
    double h_Lx = u.getLx();
    double h_Ly = u.getLy();
    double h_Lz = u.getLz();
    double h_a_t = u.a_t();

    double h_shiftX = b.getShiftX();
    double h_shiftY = b.getShiftY();
    double h_shiftZ = b.getShiftZ();

    int h_sizeI = b.getSizeI();
    int h_sizeJ = b.getSizeJ();
    int h_sizeK = b.getSizeK();

    int h_minI = b.getMinI();
    int h_minJ = b.getMinJ();
    int h_minK = b.getMinK();

    int h_size = h_sizeI * h_sizeJ * h_sizeK;

    int h_sizeBI = 1;
    int h_sizeBJ = 1;
    int h_sizeBK = h_sizeK - 2;

    double* d_block;

    SAFE_CALL(cudaMalloc((void**)&d_block, sizeof(double) * h_size));

    cudaMemcpyToSymbol(Lx, &h_Lx, sizeof(double));
    cudaMemcpyToSymbol(Ly, &h_Ly, sizeof(double));
    cudaMemcpyToSymbol(Lz, &h_Lz, sizeof(double));
    cudaMemcpyToSymbol(a_t, &h_a_t, sizeof(double));

    cudaMemcpyToSymbol(shiftX, &h_shiftX, sizeof(double));
    cudaMemcpyToSymbol(shiftY, &h_shiftY, sizeof(double));
    cudaMemcpyToSymbol(shiftZ, &h_shiftZ, sizeof(double));
    
    cudaMemcpyToSymbol(sizeI, &h_sizeI, sizeof(int));
    cudaMemcpyToSymbol(sizeJ, &h_sizeJ, sizeof(int));
    cudaMemcpyToSymbol(sizeK, &h_sizeK, sizeof(int));

    SAFE_CALL(cudaMemcpyToSymbol(minI, &h_minI, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(minJ, &h_minJ, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(minK, &h_minK, sizeof(int)));

    dim3 grid = dim3((h_sizeK - 2) / h_sizeBK, (h_sizeJ - 2) / h_sizeBJ, (h_sizeI - 2) / h_sizeBI);
    dim3 block = dim3(h_sizeBK, h_sizeBJ, h_sizeBI);

    u0<<<grid, block>>>(d_block);

    SAFE_CALL(cudaMemcpy(b.getData().data(), d_block, sizeof(double) * h_size, cudaMemcpyDeviceToHost));
}

__device__
double getElem(double* block, int i, int j, int k) {
    return block[i * (sizeJ * sizeK) + j * sizeK + k];
}

__device__
double lap_h(double* block, int i, int j, int k) {
    return (getElem(block, i - 1, j, k) - 2 * getElem(block, i, j, k) + getElem(block, i + 1, j, k)) / pow(shiftX, 2) + 
        (getElem(block, i, j - 1, k) - 2 * getElem(block, i, j, k) + getElem(block, i, j + 1, k)) / pow(shiftY, 2) +
        (getElem(block, i, j, k - 1) - 2 * getElem(block, i, j, k) + getElem(block, i, j, k + 1)) / pow(shiftZ, 2); 
}

__global__
void u1(double* block, double* u0, double tau) {
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    block[i * (sizeJ * sizeK) + j * (sizeK) + k] = u0[i * (sizeJ * sizeK) + j * (sizeK) + k] + 
        (pow(tau, 2) / 2) * lap_h(u0, i, j, k);
}

void init_u1(Block &b, const Block &u0, double tau, Function3D &u) {
    double h_shiftX = b.getShiftX();
    double h_shiftY = b.getShiftY();
    double h_shiftZ = b.getShiftZ();

    int h_sizeI = b.getSizeI();
    int h_sizeJ = b.getSizeJ();
    int h_sizeK = b.getSizeK();

    int h_size = h_sizeI * h_sizeJ * h_sizeK;

    int h_sizeBI = 1;
    int h_sizeBJ = 1;
    int h_sizeBK = h_sizeK - 2;

    dim3 grid = dim3((h_sizeK - 2) / h_sizeBK, (h_sizeJ - 2) / h_sizeBJ, (h_sizeI - 2) / h_sizeBI);
    dim3 block = dim3(h_sizeBK, h_sizeBJ, h_sizeBI);

    double* d_block;
    double* d_u0;

    SAFE_CALL(cudaMalloc((void**)&d_block, sizeof(double) * h_size));
    SAFE_CALL(cudaMalloc((void**)&d_u0, sizeof(double) * h_size));

    SAFE_CALL(cudaMemcpy(d_u0, u0.getValData().data(), sizeof(double) * h_size, cudaMemcpyHostToDevice));

    cudaMemcpyToSymbol(shiftX, &h_shiftX, sizeof(double));
    cudaMemcpyToSymbol(shiftY, &h_shiftY, sizeof(double));
    cudaMemcpyToSymbol(shiftZ, &h_shiftZ, sizeof(double));
    
    cudaMemcpyToSymbol(sizeI, &h_sizeI, sizeof(int));
    cudaMemcpyToSymbol(sizeJ, &h_sizeJ, sizeof(int));
    cudaMemcpyToSymbol(sizeK, &h_sizeK, sizeof(int));

    u1<<<grid, block>>>(d_block, d_u0, tau);

    SAFE_CALL(cudaMemcpy(b.getData().data(), d_block, sizeof(double) * h_size, cudaMemcpyDeviceToHost));
}

__global__
void global_step(double* block, double* u1, double* u0, double tau) {
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    block[i * (sizeJ * sizeK) + j * (sizeK) + k] = 2 * getElem(u1, i, j, k) - getElem(u0, i, j, k) + 
        pow(tau, 2) * lap_h(u1, i, j, k);
}

void step(Block &b, const Block& u1, const Block& u0, double tau, Function3D &u) {
    double h_shiftX = b.getShiftX();
    double h_shiftY = b.getShiftY();
    double h_shiftZ = b.getShiftZ();

    int h_sizeI = b.getSizeI();
    int h_sizeJ = b.getSizeJ();
    int h_sizeK = b.getSizeK();

    int h_size = h_sizeI * h_sizeJ * h_sizeK;

    int h_sizeBI = 1;
    int h_sizeBJ = 1;
    int h_sizeBK = h_sizeK - 2;

    dim3 grid = dim3((h_sizeK - 2) / h_sizeBK, (h_sizeJ - 2) / h_sizeBJ, (h_sizeI - 2) / h_sizeBI);
    dim3 block = dim3(h_sizeBK, h_sizeBJ, h_sizeBI);

    double* d_block;
    double* d_u1;
    double* d_u0;

    SAFE_CALL(cudaMalloc((void**)&d_block, sizeof(double) * h_size));
    SAFE_CALL(cudaMalloc((void**)&d_u1, sizeof(double) * h_size));
    SAFE_CALL(cudaMalloc((void**)&d_u0, sizeof(double) * h_size));

    SAFE_CALL(cudaMemcpy(d_u1, u1.getValData().data(), sizeof(double) * h_size, cudaMemcpyHostToDevice));
    SAFE_CALL(cudaMemcpy(d_u0, u0.getValData().data(), sizeof(double) * h_size, cudaMemcpyHostToDevice));

    cudaMemcpyToSymbol(shiftX, &h_shiftX, sizeof(double));
    cudaMemcpyToSymbol(shiftY, &h_shiftY, sizeof(double));
    cudaMemcpyToSymbol(shiftZ, &h_shiftZ, sizeof(double));
    
    cudaMemcpyToSymbol(sizeI, &h_sizeI, sizeof(int));
    cudaMemcpyToSymbol(sizeJ, &h_sizeJ, sizeof(int));
    cudaMemcpyToSymbol(sizeK, &h_sizeK, sizeof(int));

    global_step<<<grid, block>>>(d_block, d_u1, d_u0, tau);

    SAFE_CALL(cudaMemcpy(b.getData().data(), d_block, sizeof(double) * h_size, cudaMemcpyDeviceToHost));
}

__global__
void calcErrorK(double* d_error, double* block, double t) {
    extern __shared__ double sdata[];

    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    int tid = threadIdx.x;
    int ind = i * (sizeJ * sizeK) + j * sizeK + k;

    sdata[tid] = fabs(block[ind] - u((i + minI) * shiftX, (j + minJ) * shiftY, (k + minK) * shiftZ, t));
    __syncthreads();

    for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
        if (tid < s) {
            sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) d_error[blockIdx.z * gridDim.y + blockIdx.y] = sdata[0];
}

__global__
void calcErrorJ(double* d_error) {
    extern __shared__ double sdata[];

    int tid = threadIdx.x;
    int ind = blockDim.x * blockIdx.y + threadIdx.x;

    sdata[tid] = d_error[ind];
    __syncthreads();

    for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
        if (tid < s) {
            sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) d_error[blockDim.x * blockIdx.y] = sdata[0];
}

__global__
void calcErrorI(double* d_error, int size) {
    extern __shared__ double sdata[];

    int tid = threadIdx.x;
    int ind = size * threadIdx.x;

    sdata[tid] = d_error[ind];
    __syncthreads();

    for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
        if (tid < s) {
            sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) d_error[0] = sdata[0];
}

double getError(const Block &b, Function3D &u, double t) {
    double h_Lx = u.getLx();
    double h_Ly = u.getLy();
    double h_Lz = u.getLz();
    double h_a_t = u.a_t();

    double h_shiftX = b.getShiftX();
    double h_shiftY = b.getShiftY();
    double h_shiftZ = b.getShiftZ();

    int h_sizeI = b.getSizeI();
    int h_sizeJ = b.getSizeJ();
    int h_sizeK = b.getSizeK();

    int h_minI = b.getMinI();
    int h_minJ = b.getMinJ();
    int h_minK = b.getMinK();

    int h_size = h_sizeI * h_sizeJ * h_sizeK;

    double h_error;

    double* d_error;
    double* d_block;

    SAFE_CALL(cudaMalloc((void**)&d_error, sizeof(double) * (h_sizeI - 2) * (h_sizeJ - 2)));
    SAFE_CALL(cudaMalloc((void**)&d_block, sizeof(double) * h_size));

    SAFE_CALL(cudaMemcpy(d_block, b.getValData().data(), sizeof(double) * h_size, cudaMemcpyHostToDevice));

    SAFE_CALL(cudaMemcpyToSymbol(Lx, &h_Lx, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(Ly, &h_Ly, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(Lz, &h_Lz, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(a_t, &h_a_t, sizeof(double)));

    SAFE_CALL(cudaMemcpyToSymbol(shiftX, &h_shiftX, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(shiftY, &h_shiftY, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(shiftZ, &h_shiftZ, sizeof(double)));
    
    SAFE_CALL(cudaMemcpyToSymbol(sizeI, &h_sizeI, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(sizeJ, &h_sizeJ, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(sizeK, &h_sizeK, sizeof(int)));

    SAFE_CALL(cudaMemcpyToSymbol(minI, &h_minI, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(minJ, &h_minJ, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(minK, &h_minK, sizeof(int)));

    int h_sizeBI = 1;
    int h_sizeBJ = 1;
    int h_sizeBK = h_sizeK - 2;

    dim3 gridK = dim3((h_sizeK - 2) / h_sizeBK, (h_sizeJ - 2) / h_sizeBJ, (h_sizeI - 2) / h_sizeBI);
    dim3 blockK = dim3(h_sizeBK, h_sizeBJ, h_sizeBI);

    // std::cout << gridK.x << ' ' << gridK.y << ' ' << gridK.z << std::endl;
    // std::cout << blockK.x << ' ' << blockK.y << ' ' << blockK.z << std::endl;

    calcErrorK<<<gridK, blockK, blockK.x * sizeof(double)>>>(d_error, d_block, t);

    h_sizeBI = 1;
    h_sizeBJ = h_sizeJ - 2;
    h_sizeBK = 1;

    dim3 gridJ = dim3((h_sizeJ - 2) / h_sizeBJ, (h_sizeI - 2) / h_sizeBI, 1);
    dim3 blockJ = dim3(h_sizeBJ, h_sizeBI, 1);

    // std::cout << gridJ.x << ' ' << gridJ.y << ' ' << gridJ.z << std::endl;
    // std::cout << blockJ.x << ' ' << blockJ.y << ' ' << blockJ.z << std::endl;

    calcErrorJ<<<gridJ, blockJ, blockJ.x * sizeof(double)>>>(d_error);

    h_sizeBI = h_sizeI - 2;
    h_sizeBJ = 1;
    h_sizeBK = 1;

    dim3 gridI = dim3((h_sizeI - 2) / h_sizeBI, 1, 1);
    dim3 blockI = dim3(h_sizeBI, 1, 1);

    calcErrorI<<<gridI, blockI, blockI.x * sizeof(double)>>>(d_error, blockJ.x);

    SAFE_CALL(cudaMemcpy(&h_error, d_error, sizeof(double), cudaMemcpyDeviceToHost));
    cudaFree(d_error);
    cudaFree(d_block);
    return h_error;
}