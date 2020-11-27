#pragma onece
#include <mpi.h>
#include <vector>

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

public:  
    Process(int countI, int countJ, int countK, int N);

    void sendDownI(std::vector<double> downI, MPI_Request &request);
    void sendUpI(std::vector<double> upI, MPI_Request &request);
    void sendDownJ(std::vector<double> downJ, MPI_Request &request);
    void sendUpJ(std::vector<double> upJ, MPI_Request &request);
    void sendDownK(std::vector<double> downK, MPI_Request &request);
    void sendUpK(std::vector<double> upK, MPI_Request &request);

    std::vector<double> recvDownI(MPI_Status &status);
    std::vector<double> recvUpI(MPI_Status &status);
    std::vector<double> recvDownJ(MPI_Status &status);
    std::vector<double> recvUpJ(MPI_Status &status);
    std::vector<double> recvDownK(MPI_Status &status);
    std::vector<double> recvUpK(MPI_Status &status);

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