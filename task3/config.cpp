#include <cstdlib>
#include <cstring>
#include <iostream>


//#include <cuda_runtime.h>

#include "config.h"

Config::Config(int argc, char** argv) 
    : _argc(argc), _argv(argv),
        _Lx(1.0f), _Ly(1.0f), _Lz(1.0f),
        _T(0.1f), _K(20), _N(128)
{
    for (int i = 1; i < argc; i++) {
        if (!strcmp("-L", argv[i])) {
            _Lx, _Ly, _Lz = atof(argv[++i]);
        } else if (!strcmp("-Lx", argv[i])) {
            _Lx = atof(argv[++i]);
        } else if (!strcmp("-Ly", argv[i])) {
            _Ly = atof(argv[++i]);
        } else if (!strcmp("-Lz", argv[i])) {
            _Lz = atof(argv[++i]);
        } else if (!strcmp("-T", argv[i])) {
            _T = atof(argv[++i]);
        } else if (!strcmp("-K", argv[i])) {
            _K = atoi(argv[++i]);
        } else if (!strcmp("-N", argv[i])) {
            _N = atoi(argv[++i]);
        }
    }
    _tau = _T / _K;
}

void Config::printConfig() {
    std::cout << "=============================================" << std::endl;
    std::cout <<  "C++:" << std::endl;
    std::cout << "\tStandart:" << __cplusplus << std::endl;
    
    // std::cout << "Cuda:" << std::endl;
    // cudaDeviceProp devProp;
    // cudaGetDeviceProperties ( &devProp, 0 );
    // std::cout << "\tCompute capability     : " << devProp.major << "." << devProp.minor << std::endl;
    // std::cout << "\tName                   : " << devProp.name << std::endl;
    // std::cout << "\tTotal Global Memory    : " << devProp.totalGlobalMem << std::endl;
    // std::cout << "\tShared memory per block: " << devProp.sharedMemPerBlock << std::endl;
    // std::cout << "\tRegisters per block    : " << devProp.regsPerBlock << std::endl;
    // std::cout << "\tWarp size              : " << devProp.warpSize << std::endl;
    // std::cout << "\tMax threads per block  : " << devProp.maxThreadsPerBlock << std::endl;
    // std::cout << "\tTotal constant memory  : " << devProp.totalConstMem << std::endl;

    std::cout << "Config:" << std::endl;
    std::cout << "\tLx: " << _Lx << std::endl;
    std::cout << "\tLy: " << _Ly << std::endl;
    std::cout << "\tLz: " << _Lz << std::endl;
    std::cout << "\tT: " << _T << std::endl;
    std::cout << "\tK: " << _K << std::endl;
    std::cout << "\ttau: " << _tau << std::endl;
    std::cout << "\tN: " << _N << std::endl;
    std::cout << "Process:" << std::endl;
    std::cout << "\tndim i: " << _ndims[0] << std::endl;
    std::cout << "\tndim j: " << _ndims[1] << std::endl;
    std::cout << "\tndim k: " << _ndims[2] << std::endl;
    std::cout << "=============================================" << std::endl;
}