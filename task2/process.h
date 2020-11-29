#pragma once
#include <mpi.h>
#include <vector>

#include "block.h"
#include "function.h"

class Process {
    int _rank;
    int _p_count;

    int _countI;
    int _countJ;
    int _countK;

    int _curI;
    int _curJ;
    int _curK;

    int _bsize_I;
    int _bsize_J;
    int _bsize_K;

    int _bmin_i;
    int _bmin_j;
    int _bmin_k;

    double _startTime;

    MPI_Request _request;
    MPI_Request _downI;
    MPI_Request _upI;
    MPI_Request _downJ;
    MPI_Request _upJ;
    MPI_Request _downK;
    MPI_Request _upK;
    MPI_Status _status;

public:  
    Process(int countI, int countJ, int countK, int N, int argc, char** argv);
    ~Process();

    void printError(const Block& b, Function3D &u, double t);

    void update(Block &b);

    void sendDownI(std::vector<double> &downI);
    void sendUpI(std::vector<double> &upI);
    void sendDownJ(std::vector<double> &downJ);
    void sendUpJ(std::vector<double> &upJ);
    void sendDownK(std::vector<double> &downK);
    void sendUpK(std::vector<double> &upK);

    std::vector<double> recvDownI();
    std::vector<double> recvUpI();
    std::vector<double> recvDownJ();
    std::vector<double> recvUpJ();
    std::vector<double> recvDownK();
    std::vector<double> recvUpK();

    int getRank() { return _rank; }
    int getOtherRank(int i, int j, int k) { return i * (_countJ * _countK) + j * _countK + k; }

    int getI() { return _curI; }
    int getJ() { return _curJ; }
    int getK() { return _curK; }

    int getSizeX() { return _bsize_I; }
    int getSizeY() { return _bsize_J; }
    int getSizeZ() { return _bsize_K; }

    int get_i() { return _bmin_i; }
    int get_j() { return _bmin_j; }
    int get_k() { return _bmin_k; }
};