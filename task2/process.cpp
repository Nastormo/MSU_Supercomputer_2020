#include "process.h"

Process::Process(int countI, int countJ, int countK, int N) 
        : _countI(countI),
        _countJ(countJ),
        _countK(countK) {
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_p_count);

    _curI = _rank / (_countJ * _countK);;
    _curJ = (_rank - (_countJ * _countK) * _curI) / _countK;
    _curK = _rank % _countK;

    _bsize_I = N / _countI;
    _bsize_J = N / _countJ;
    _bsize_K = N / _countK;

    _bmin_i = _bsize_I * _curI;
    _bmin_j = _bsize_J * _curJ;
    _bmin_k = _bsize_K * _curK;
}

void Process::sendDownI(std::vector<double> downI, MPI_Request &request) {
    if (_curI > 0) {
        MPI_Isend(downI.data(), downI.size(), MPI_DOUBLE, 
            getOtherRank(_curI - 1, _curJ, _curK), 
            0, MPI_COMM_WORLD, &request);
    }
}

void Process::sendUpI(std::vector<double> upI, MPI_Request &request) {
    if (_curI + 1 < _countI) {
        MPI_Isend(upI.data(), upI.size(), MPI_DOUBLE, 
            getOtherRank(_curI + 1, _curJ, _curK), 
            0, MPI_COMM_WORLD, &request);
    }
}

void Process::sendDownJ(std::vector<double> downJ, MPI_Request &request) {
    if (_curJ > 0) {
        MPI_Isend(downJ.data(), downJ.size(), MPI_DOUBLE, 
            getOtherRank(_curI, _curJ - 1, _curK), 
            0, MPI_COMM_WORLD, &request);
    }
}

void Process::sendUpJ(std::vector<double> upJ, MPI_Request &request) {
    if (_curJ + 1 < _countJ) {
        MPI_Isend(upJ.data(), upJ.size(), MPI_DOUBLE, 
            getOtherRank(_curI, _curJ + 1, _curK), 
            0, MPI_COMM_WORLD, &request);
    }
}

void Process::sendDownK(std::vector<double> downK, MPI_Request &request) {
    if (_curK > 0) {
        MPI_Isend(downK.data(), downK.size(), MPI_DOUBLE, 
            getOtherRank(_curI, _curJ, _curK - 1), 
            0, MPI_COMM_WORLD, &request);
    }
}

void Process::sendUpK(std::vector<double> upK, MPI_Request &request) {
    if (_curK + 1 < _countK) {
        MPI_Isend(upK.data(), upK.size(), MPI_DOUBLE, 
            getOtherRank(_curI, _curJ, _curK + 1), 
            0, MPI_COMM_WORLD, &request);
    }
}

//Recv
std::vector<double> Process::recvDownI(MPI_Status &status) {
    std::vector<double> downI (_bsize_I, 0);
    if (_curI > 0) {
        MPI_Recv(downI.data(), downI.size(), MPI_DOUBLE,
            getOtherRank(_curI - 1, _curJ, _curK),
            0, MPI_COMM_WORLD, &status);
    }
    return downI;
}

std::vector<double> Process::recvUpI(MPI_Status &status) {
    std::vector<double> upI (_bsize_I, 0);
    if (_curI + 1 < _countI) {
        MPI_Recv(upI.data(), upI.size(), MPI_DOUBLE,
            getOtherRank(_curI + 1, _curJ, _curK),
            0, MPI_COMM_WORLD, &status);
    }
    return upI;
}

std::vector<double> Process::recvDownJ(MPI_Status &status) {
    std::vector<double> downJ (_bsize_J, 0);
    if (_curJ > 0) {
        MPI_Recv(downJ.data(), downJ.size(), MPI_DOUBLE,
            getOtherRank(_curI, _curJ - 1, _curK),
            0, MPI_COMM_WORLD, &status);
    }
    return downJ;
}

std::vector<double> Process::recvUpJ(MPI_Status &status) {
    std::vector<double> upJ (_bsize_J, 0);
    if (_curJ + 1 < _countJ) {
        MPI_Recv(upJ.data(), upJ.size(), MPI_DOUBLE,
            getOtherRank(_curI, _curJ + 1, _curK),
            0, MPI_COMM_WORLD, &status);
    }
    return upJ;
}

std::vector<double> Process::recvDownK(MPI_Status &status) {
    std::vector<double> downK (_bsize_K, 0);
    if (_curK > 0) {
        MPI_Recv(downK.data(), downK.size(), MPI_DOUBLE,
            getOtherRank(_curI, _curJ, _curK - 1),
            0, MPI_COMM_WORLD, &status);
    }
    return downK;
}

std::vector<double> Process::recvUpK(MPI_Status &status) {
    std::vector<double> upK (_bsize_K, 0);
    if (_curK + 1 < _countK) {
        MPI_Recv(upK.data(), upK.size(), MPI_DOUBLE,
            getOtherRank(_curI, _curJ, _curK + 1),
            0, MPI_COMM_WORLD, &status);
    }
    return upK;
}