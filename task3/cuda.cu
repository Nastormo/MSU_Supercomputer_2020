#include "parallel.h"

#include "cuda_common.h"
#include <cuda.h>

#include <vector>
#include <iostream>
#include <math.h>
#include <mpi.h>

__constant__ double a_t;
__constant__ double L[3];
__constant__ double bshift[3];
__constant__ int bsize[3];
__constant__ int bmin[3];

__device__
double d_cu(double x, double y, double z, double t) {
    return sin((M_PI / L[0]) * x) * 
           sin((M_PI / L[1]) * y) * 
           sin((M_PI / L[2]) * z) *
           cos(a_t * t);;
}

__device__
double d_getElem(double* block, int i, int j, int k) {
    return block[i * (bsize[1] * bsize[2]) + j * bsize[2] + k];
}

__device__
double d_lap_h(double* block, int i, int j, int k) {
    double central = d_getElem(block, i, j, k);
    return (d_getElem(block, i - 1, j, k) - 2 * central + d_getElem(block, i + 1, j, k)) / pow(bshift[0], 2) + 
           (d_getElem(block, i, j - 1, k) - 2 * central + d_getElem(block, i, j + 1, k)) / pow(bshift[1], 2) +
           (d_getElem(block, i, j, k - 1) - 2 * central + d_getElem(block, i, j, k + 1)) / pow(bshift[2], 2); 
}

__global__
void g_u0(double* block)
{
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    int ind = i * (bsize[1] * bsize[2]) + j * bsize[2] + k;
    block[ind] = 
        d_cu((i + bmin[0]) * bshift[0], (j + bmin[1]) * bshift[1], (k + bmin[2]) * bshift[2], 0);
}

__global__
void g_u1(double* block, double* u0, double tau) {
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    int ind = i * (bsize[1] * bsize[2]) + j * bsize[2] + k;
    block[ind] = u0[ind] + (pow(tau, 2) / 2) * d_lap_h(u0, i, j, k);
}

__global__
void g_step(double* block, double* u1, double* u0, double tau) {
    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    int ind = i * (bsize[1] * bsize[2]) + j * bsize[2] + k;
    block[ind] = 2 * u1[ind] - u0[ind] + pow(tau, 2) * d_lap_h(u1, i, j, k);
}

Parallel::Parallel(Config& conf, Function3D& u)
    : _conf(conf), _p(conf), _u(u), _t(0), _tau(conf.getTau()), _K(conf.getK())
{
    std::vector<double>L = _u.getL();
    int N = _conf.getN();

    std::vector<int> bsize = _p.getSize();
    std::vector<int> bmin = _p.getMin();
    std::vector<double> bshift = { L[0]/(N + 1), L[1]/(N + 1), L[2]/(N + 1) };
    
    for (int i = 0; i < 3; i++) {
        _massB.push_back(Block(bsize, bmin, bshift));
        _massB[i].init_cuda();
    }
    _bshift = _massB[0].getShift();
    _bsize = _massB[0].getSize();
    _bmin = _massB[0].getMin();
}

void Parallel::init_u0() {
    double h_a_t = _u.a_t();
    std::vector<double> h_L = _u.getL();

    int h_size = _bsize[0] * _bsize[1] * _bsize[2];

    double* d_u0 = _massB[_t % 3].getDataLink();

    checkCudaErrors();
    SAFE_CALL(cudaMemcpyToSymbol(a_t, &h_a_t, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(L, h_L.data(), sizeof(double) * 3));
    SAFE_CALL(cudaMemcpyToSymbol(bshift, _bshift.data(), sizeof(double) * 3));
    SAFE_CALL(cudaMemcpyToSymbol(bsize, _bsize.data(), sizeof(int) * 3));
    SAFE_CALL(cudaMemcpyToSymbol(bmin, _bmin.data(), sizeof(int) * 3));

    dim3 grid = dim3(1, (_bsize[1] - 2), (_bsize[0] - 2));
    dim3 block = dim3(_bsize[2] - 2, 1, 1);

    g_u0<<<grid, block>>>(d_u0);

    checkCudaErrors();
}

void Parallel::init_u1() {
    double* d_u0 = _massB[(_t - 1) % 3].getDataLink();
    double* d_u1 = _massB[_t % 3].getDataLink();

    cudaMemcpyToSymbol(bshift, _bshift.data(), sizeof(double) * 3);
    cudaMemcpyToSymbol(bsize, _bsize.data(), sizeof(int) * 3);

    dim3 grid = dim3(1, (_bsize[1] - 2), (_bsize[0] - 2));
    dim3 block = dim3(_bsize[2] - 2, 1, 1);

    g_u1<<<grid, block>>>(d_u1, d_u0, _tau);

    checkCudaErrors();
}

void Parallel::step() {
    double* d_u0 = _massB[(_t-2) % 3].getDataLink();
    double* d_u1 = _massB[(_t-1) % 3].getDataLink();
    double* d_u2 = _massB[(_t) % 3].getDataLink();

    cudaMemcpyToSymbol(bshift, _bshift.data(), sizeof(double) * 3);
    cudaMemcpyToSymbol(bsize, _bsize.data(), sizeof(int) * 3);

    dim3 grid = dim3(1, (_bsize[1] - 2), (_bsize[0] - 2));
    dim3 block = dim3(_bsize[2] - 2, 1, 1);

    g_step<<<grid, block>>>(d_u2, d_u1, d_u0, _tau);

    checkCudaErrors();
}

void Parallel::printError() {
    double time = _t * _tau;
    _p.printError(time, _massB[_t % 3].getError(_u, time));
}

void Parallel::update() {
    _p.update(_massB[_t % 3]);
}

void Parallel::process() {
    init_u0();
    printError();
    update();
    _t++;
    
    init_u1();
    printError();
    update();
    _t++;

    while (_t < _K) {
        step();
        printError();
        update();
        _t++;
    }
}