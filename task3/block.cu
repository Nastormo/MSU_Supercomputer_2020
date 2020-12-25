#include "block.h"
#include "cuda_common.h"
#include <cuda.h>

__constant__ double a_t;
__constant__ double L[3];
__constant__ double bshift[3];
__constant__ int bsize[3];
__constant__ int bmin[3];

__device__
double d_bu(double x, double y, double z, double t) {
    return sin((M_PI / L[0]) * x) * 
           sin((M_PI / L[1]) * y) * 
           sin((M_PI / L[2]) * z) *
           cos(a_t * t);;
}

__global__
void g_calcErrorK(double* d_error, double* block, double t) {
    extern __shared__ double sdata[];

    int i = blockIdx.z + 1;
    int j = blockIdx.y + 1;
    int k = threadIdx.x + 1;
    int tid = threadIdx.x;
    int ind = i * (bsize[1] * bsize[2]) + j * bsize[2] + k;

    sdata[tid] = fabs(block[ind] - d_bu((i + bmin[0]) * bshift[0], (j + bmin[1]) * bshift[1], (k + bmin[2]) * bshift[2], t));
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
void g_calcErrorJ(double* d_errorJ, double* d_errorK) {
    extern __shared__ double sdata[];

    int tid = threadIdx.x;
    int ind = blockDim.x * blockIdx.y + threadIdx.x;

    sdata[tid] = d_errorK[ind];
    __syncthreads();

    for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
        if (tid < s) {
            sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) d_errorJ[blockIdx.y] = sdata[0];
}

__global__
void g_calcErrorI(double* d_errorJ) {
    extern __shared__ double sdata[];

    int tid = threadIdx.x;

    sdata[tid] = d_errorJ[tid];
    __syncthreads();

    for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
        if (tid < s) {
            sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) d_errorJ[0] = sdata[0];
}

__global__
void g_getSlice(double* block, double* slice, int axis, int item) {
    int ijk[3];
    ijk[(axis + 2) % 3] = threadIdx.x + 1;
    ijk[(axis + 1) % 3] = blockIdx.y + 1;
    ijk[axis] = item;

    slice[threadIdx.x * gridDim.y + blockIdx.y] = 
        block[ijk[0] * (bsize[1] * bsize[2]) + ijk[1] * bsize[2] + ijk[2]];
}

__global__
void g_setSlice(double* block, double* slice, int axis, int item) {
    int ijk[3];
    ijk[(axis + 2) % 3] = threadIdx.x + 1;
    ijk[(axis + 1) % 3] = blockIdx.y + 1;
    ijk[axis] = item;

    block[ijk[0] * (bsize[1] * bsize[2]) + ijk[1] * bsize[2] + ijk[2]] = 
        slice[threadIdx.x * gridDim.y + blockIdx.y];     
}

Block::Block(std::vector<int>& size, std::vector<int>& min, 
    std::vector<double>& shift) 
        : _size({size[0] + 2, size[1] + 2, size[2] + 2}), 
          _min(min), _shift(shift)
{
    _size = {size[0] + 2, size[1] + 2, size[2] + 2};
    _shift = shift;
    _raw.resize(_size[0] * _size[1] * _size[2], 0.0);
}

void Block::init_cuda() {
    int max_ssize = std::max((_size[0] - 2) * (_size[1] - 2), 
        std::max((_size[0] - 2) * (_size[2] - 2), (_size[1] - 2) * (_size[2] - 2)));
    SAFE_CALL(cudaMalloc((void**)&_d_raw, sizeof(double) * _size[0] * _size[1] * _size[2]));
    SAFE_CALL(cudaMalloc((void**)&_d_errorK, sizeof(double) * (_size[0] - 2) * (_size[1] - 2)));
    SAFE_CALL(cudaMalloc((void**)&_d_errorJ, sizeof(double) * (_size[0] - 2)));
    SAFE_CALL(cudaMalloc((void**)&_d_slice, sizeof(double) * max_ssize));
}

void Block::destroy_cuda() {
    cudaFree(_d_raw);
    cudaFree(_d_errorK);
    cudaFree(_d_errorJ);
    cudaFree(_d_slice);
}

////////////////////////////////////////////////////////////
//TODO
void Block::printBlock() const { return; };
void Block::printDiff(Function3D &u, double t) const { return; };
void Block::saveBlock(std::string &str) const { return; };
////////////////////////////////////////////////////////////

double Block::getError(Function3D &u, double t) const
{
    double h_a_t = u.a_t();
    std::vector<double> h_L = u.getL();;

    double h_error;

    dim3 grid;
    dim3 block;

    checkCudaErrors();
    SAFE_CALL(cudaMemcpyToSymbol(a_t, &h_a_t, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(L, h_L.data(), sizeof(double) * 3));
    SAFE_CALL(cudaMemcpyToSymbol(bshift, _shift.data(), sizeof(double) * 3));
    SAFE_CALL(cudaMemcpyToSymbol(bsize, _size.data(), sizeof(int) * 3));
    SAFE_CALL(cudaMemcpyToSymbol(bmin, _min.data(), sizeof(int) * 3));

    grid = dim3(1, (_size[1] - 2), (_size[0] - 2));
    block = dim3(_size[2] - 2, 1, 1);

    g_calcErrorK<<<grid, block, block.x * sizeof(double)>>>(_d_errorK, _d_raw, t);
    checkCudaErrors();

    grid = dim3(1, (_size[0] - 2), 1);
    block = dim3(_size[1] - 2, 1, 1);

    g_calcErrorJ<<<grid, block, block.x * sizeof(double)>>>(_d_errorJ, _d_errorK);
    checkCudaErrors();

    grid = dim3(1, 1, 1);
    block = dim3(_size[0] - 2, 1, 1);

    g_calcErrorI<<<grid, block, block.x * sizeof(double)>>>(_d_errorJ);
    checkCudaErrors();

    SAFE_CALL(cudaMemcpy(&h_error, _d_errorJ, sizeof(double), cudaMemcpyDeviceToHost));
    return h_error;
}

double& Block::getElem(int i, int j, int k) {
    return _raw[i * (_size[1] * _size[2]) + j * _size[2] + k];
}

double Block::getValElem(int i, int j, int k) const {
    return _raw[i * (_size[1] * _size[2]) + j * _size[2] + k];
}

std::vector<double> Block::getSlice(int axis, int item) const {
    int ssize = (_size[(axis + 1) % 3] - 2) * (_size[(axis + 2) % 3] - 2);

    SAFE_CALL(cudaMemcpyToSymbol(bsize, _size.data(), sizeof(int) * 3));

    dim3 grid = dim3(1, _size[(axis + 1) % 3] - 2, 1);
    dim3 block = dim3(_size[(axis + 2) % 3] - 2, 1, 1);

    g_getSlice<<<grid, block>>>(_d_raw, _d_slice, axis, item);

    std::vector<double> h_slice (ssize, 0);

    SAFE_CALL(cudaMemcpy(h_slice.data(), _d_slice, sizeof(double) * ssize, cudaMemcpyDeviceToHost));

    return h_slice;
};

void Block::setSlice(const std::vector<double>& slice, int axis, int item) {
    int ssize = slice.size();

    SAFE_CALL(cudaMemcpyToSymbol(bsize, _size.data(), sizeof(int) * 3));
    SAFE_CALL(cudaMemcpy(_d_slice, slice.data(), sizeof(double) * ssize, cudaMemcpyHostToDevice));

    dim3 grid = dim3(1, _size[(axis + 1) % 3] - 2, 1);
    dim3 block = dim3(_size[(axis + 2) % 3] - 2, 1, 1);

    g_setSlice<<<grid, block>>>(_d_raw, _d_slice, axis, item);
}

double Block::lap_h(int i, int j, int k) const {
    return (getValElem(i - 1, j, k) - 2 * getValElem(i, j, k) + getValElem(i + 1, j, k)) / pow(_shift[0], 2) + 
        (getValElem(i, j - 1, k) - 2 * getValElem(i, j, k) + getValElem(i, j + 1, k)) / pow(_shift[1], 2) +
        (getValElem(i, j, k - 1) - 2 * getValElem(i, j, k) + getValElem(i, j, k + 1)) / pow(_shift[2], 2); 
}