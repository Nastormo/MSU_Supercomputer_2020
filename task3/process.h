#pragma once
#include <mpi.h>
#include <vector>
#include <iostream>

#include "config.h"
#include "block.h"
#include "function.h"


class Process {
    int _rank;
    int _p_count;

    std::vector<int> _ndims;
    std::vector<int> _position;
    std::vector<int> _bsize;
    std::vector<int> _bmin;

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
    Process(Config conf);
    ~Process();

    void printError(double t, double otherError);

    void update(Block &b);

    void Send(std::vector<double>& slice, int axis, int shift);
    std::vector<double> Recv(int axis, int shift);

    int getRank() { return _rank; }
    int getPCount() { return _p_count; } 
    int getOtherRank(int i, int j, int k) { return i * (_ndims[1] * _ndims[2]) + j * _ndims[2] + k; }
    int getOtherRank(const std::vector<int> &position);

    int getI() { return _position[0]; }
    int getJ() { return _position[1]; }
    int getK() { return _position[2]; }

    std::vector<int> getSize() { return _bsize; }
    int getSizeI() { return _bsize[0]; }
    int getSizeJ() { return _bsize[1]; }
    int getSizeK() { return _bsize[2]; }

    std::vector<int> getMin() { return _bmin; }
    int get_i() { return _bmin[0]; }
    int get_j() { return _bmin[1]; }
    int get_k() { return _bmin[2]; }
};