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

__constant__ double a_t;
__constant__ double L[3];
__constant__ double bshift[3];
__constant__ int bsize[3];
__constant__ int bmin[3];

__device__
double u(double x, double y, double z, double t) {
    return sin((M_PI / L[0]) * x) * 
           sin((M_PI / L[1]) * y) * 
           sin((M_PI / L[2]) * z) *
           cos(a_t * t);;
}

__global__
void u0(double* block)
{
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    block[i * (bsize[1] * bsize[2]) + j * (bsize[2]) + k] = 
        u((i + bmin[0]) * bshift[0], (j + bmin[1]) * bshift[1], (k + bmin[2]) * bshift[2], 0.0f);
}

void init_u0(Block &b, Function3D &u) {
    double h_a_t = u.a_t();
    std::vector<double> h_L = u.getL();
    std::vector<double> h_bshift = b.getShift();
    std::vector<int> h_bsize = b.getSize();
    std::vector<int> h_bmin = b.getMin();

    int h_size = h_bsize[0] * h_bsize[1] * h_bsize[2];

    double* d_block;

    SAFE_CALL(cudaMalloc((void**)&d_block, sizeof(double) * h_size));

    cudaMemcpyToSymbol(a_t, &h_a_t, sizeof(double));
    cudaMemcpyToSymbol(L, h_L.data(), sizeof(double) * 3);
    cudaMemcpyToSymbol(bshift, h_bshift.data(), sizeof(double) * 3);
    cudaMemcpyToSymbol(bsize, h_bsize.data(), sizeof(int) * 3);
    SAFE_CALL(cudaMemcpyToSymbol(bmin, h_bmin.data(), sizeof(int) * 3));

    dim3 grid = dim3(1, (h_bsize[1] - 2), (h_bsize[0] - 2));
    dim3 block = dim3(h_bsize[2] - 2, 1, 1);

    u0<<<grid, block>>>(d_block);

    SAFE_CALL(cudaMemcpy(b.getData().data(), d_block, sizeof(double) * h_size, cudaMemcpyDeviceToHost));
    cudaFree(d_block);
}

__device__
double getElem(double* block, int i, int j, int k) {
    return block[i * (bsize[1] * bsize[2]) + j * bsize[2] + k];
}

__device__
double lap_h(double* block, int i, int j, int k) {
    double central = getElem(block, i, j, k);
    return (getElem(block, i - 1, j, k) - 2 * central + getElem(block, i + 1, j, k)) / pow(bshift[0], 2) + 
           (getElem(block, i, j - 1, k) - 2 * central + getElem(block, i, j + 1, k)) / pow(bshift[1], 2) +
           (getElem(block, i, j, k - 1) - 2 * central + getElem(block, i, j, k + 1)) / pow(bshift[2], 2); 
}

__global__
void u1(double* block, double* u0, double tau) {
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    int ind = i * (bsize[1] * bsize[2]) + j * bsize[2] + k;
    block[ind] = u0[ind] + (pow(tau, 2) / 2) * lap_h(u0, i, j, k);
}

void init_u1(Block &b, const Block &u0, double tau, Function3D &u) {
    std::vector<double> h_bshift = b.getShift();
    std::vector<int> h_bsize = b.getSize();

    int h_size = h_bsize[0] * h_bsize[1] * h_bsize[2];

    double* d_block;
    double* d_u0;

    SAFE_CALL(cudaMalloc((void**)&d_block, sizeof(double) * h_size));
    SAFE_CALL(cudaMalloc((void**)&d_u0, sizeof(double) * h_size));

    SAFE_CALL(cudaMemcpy(d_u0, u0.getValData().data(), sizeof(double) * h_size, cudaMemcpyHostToDevice));

    cudaMemcpyToSymbol(bshift, h_bshift.data(), sizeof(double) * 3);
    cudaMemcpyToSymbol(bsize, h_bsize.data(), sizeof(int) * 3);

    dim3 grid = dim3(1, (h_bsize[1] - 2), (h_bsize[0] - 2));
    dim3 block = dim3(h_bsize[2] - 2, 1, 1);

    u1<<<grid, block>>>(d_block, d_u0, tau);

    SAFE_CALL(cudaMemcpy(b.getData().data(), d_block, sizeof(double) * h_size, cudaMemcpyDeviceToHost));
    cudaFree(d_block);
    cudaFree(d_u0);
}

__global__
void global_step(double* block, double* u1, double* u0, double tau) {
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    int ind = i * (bsize[1] * bsize[2]) + j * bsize[2] + k;
    block[ind] = 2 * u1[ind] - u0[ind] + pow(tau, 2) * lap_h(u1, i, j, k);
}

void step(Block &b, const Block& u1, const Block& u0, double tau, Function3D &u) {
    std::vector<double> h_bshift = b.getShift();
    std::vector<int> h_bsize = b.getSize();

    int h_size = h_bsize[0] * h_bsize[1] * h_bsize[2];

    double* d_block;
    double* d_u1;
    double* d_u0;

    SAFE_CALL(cudaMalloc((void**)&d_block, sizeof(double) * h_size));
    SAFE_CALL(cudaMalloc((void**)&d_u1, sizeof(double) * h_size));
    SAFE_CALL(cudaMalloc((void**)&d_u0, sizeof(double) * h_size));

    SAFE_CALL(cudaMemcpy(d_u1, u1.getValData().data(), sizeof(double) * h_size, cudaMemcpyHostToDevice));
    SAFE_CALL(cudaMemcpy(d_u0, u0.getValData().data(), sizeof(double) * h_size, cudaMemcpyHostToDevice));

    cudaMemcpyToSymbol(bshift, h_bshift.data(), sizeof(double) * 3);
    cudaMemcpyToSymbol(bsize, h_bsize.data(), sizeof(int) * 3);

    dim3 grid = dim3(1, (h_bsize[1] - 2), (h_bsize[0] - 2));
    dim3 block = dim3(h_bsize[2] - 2, 1, 1);

    global_step<<<grid, block>>>(d_block, d_u1, d_u0, tau);

    SAFE_CALL(cudaMemcpy(b.getData().data(), d_block, sizeof(double) * h_size, cudaMemcpyDeviceToHost));
    cudaFree(d_block);
    cudaFree(d_u1);
    cudaFree(d_u0);
}

__global__
void calcErrorK(double* d_error, double* block, double t) {
    extern __shared__ double sdata[];

    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    int tid = threadIdx.x;
    int ind = i * (bsize[1] * bsize[2]) + j * bsize[2] + k;

    sdata[tid] = fabs(block[ind] - u((i + bmin[0]) * bshift[0], (j + bmin[1]) * bshift[1], (k + bmin[2]) * bshift[2], t));
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
    double h_a_t = u.a_t();
    std::vector<double> h_L = u.getL();
    std::vector<double> h_bshift = b.getShift();
    std::vector<int> h_bsize = b.getSize();
    std::vector<int> h_bmin = b.getMin();

    int h_size = h_bsize[0] * h_bsize[1] * h_bsize[2];

    double h_error;

    dim3 grid;
    dim3 block;

    double* d_error;
    double* d_block;

    SAFE_CALL(cudaMalloc((void**)&d_error, sizeof(double) * (h_bsize[0] - 2) * (h_bsize[1] - 2)));
    SAFE_CALL(cudaMalloc((void**)&d_block, sizeof(double) * h_size));

    SAFE_CALL(cudaMemcpy(d_block, b.getValData().data(), sizeof(double) * h_size, cudaMemcpyHostToDevice));

    cudaMemcpyToSymbol(a_t, &h_a_t, sizeof(double));
    cudaMemcpyToSymbol(L, h_L.data(), sizeof(double) * 3);
    cudaMemcpyToSymbol(bshift, h_bshift.data(), sizeof(double) * 3);
    cudaMemcpyToSymbol(bsize, h_bsize.data(), sizeof(int) * 3);
    SAFE_CALL(cudaMemcpyToSymbol(bmin, h_bmin.data(), sizeof(int) * 3));

    grid = dim3(1, (h_bsize[1] - 2), (h_bsize[0] - 2));
    block = dim3(h_bsize[2] - 2, 1, 1);

    calcErrorK<<<grid, block, block.x * sizeof(double)>>>(d_error, d_block, t);

    grid = dim3(1, (h_bsize[0] - 2), 1);
    block = dim3(h_bsize[1] - 2, 1, 1);

    calcErrorJ<<<grid, block, block.x * sizeof(double)>>>(d_error);

    grid = dim3(1, 1, 1);
    block = dim3(h_bsize[0] - 2, 1, 1);

    calcErrorI<<<grid, block, block.x * sizeof(double)>>>(d_error, h_bsize[1] - 2);

    SAFE_CALL(cudaMemcpy(&h_error, d_error, sizeof(double), cudaMemcpyDeviceToHost));
    cudaFree(d_error);
    cudaFree(d_block);
    return h_error;
}