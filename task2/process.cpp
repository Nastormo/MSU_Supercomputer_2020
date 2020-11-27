#include "process.h"

Process::Process(int countI, int countJ, int countK, int N, MPI_Request &request, MPI_Status &status) 
        : _countI(countI),
        _countJ(countJ),
        _countK(countK),
        _request(request),
        _status(status) {
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_p_count);

    _curI = _rank / (_countJ * _countK);;
    _curJ = (_rank - (_countJ * _countK) * _curI) / _countK;
    _curK = _rank % _countK;

    _bsize_I = (N - 1) / _countI;
    _bsize_J = (N - 1) / _countJ;
    _bsize_K = (N - 1) / _countK;

    _bmin_i = _bsize_I * _curI;
    _bmin_j = _bsize_J * _curJ;
    _bmin_k = _bsize_K * _curK;
}

void Process::printError(const Block& b, Function3D &u, double t) {
    double error = b.getError(u, t);
    if (_rank == 0) {
        double otherError;
        for (int i = 1; i < _p_count; i++) {
            MPI_Recv(&otherError, 1, MPI_DOUBLE,
                i, 1, MPI_COMM_WORLD, &_status);
            error = std::max(otherError, error);
        }
        printf("t: %lf Error: %lf\n", t, error);
    } else {
        //printf("Other t: %lf Error: %lf\n", t, error);
        MPI_Isend(&error, 1, MPI_DOUBLE, 
            0, 1, MPI_COMM_WORLD, &_request);
    }
}

void Process::update(Block &b) {
    sendNeighbors(b);
    recvNeighbors(b);
}

void Process::sendNeighbors(const Block &b) {
    sendDownI(b.getDownI());
    sendUpI(b.getUpI());
    sendDownJ(b.getDownJ());
    sendUpJ(b.getUpJ());
    sendDownK(b.getDownK());
    sendUpK(b.getUpK());
}

void Process::recvNeighbors(Block &b) {
    b.setDownI(recvDownI());
    b.setUpI(recvUpI());
    b.setDownJ(recvDownJ());
    b.setUpJ(recvUpJ());
    b.setDownK(recvDownK());
    b.setUpK(recvUpK());
}

void Process::sendDownI(std::vector<double> downI) {
    if (_curI > 0) {
        MPI_Isend(downI.data(), downI.size(), MPI_DOUBLE, 
            getOtherRank(_curI - 1, _curJ, _curK), 
            0, MPI_COMM_WORLD, &_request);
    }
}

void Process::sendUpI(std::vector<double> upI) {
    if (_curI + 1 < _countI) {
        MPI_Isend(upI.data(), upI.size(), MPI_DOUBLE, 
            getOtherRank(_curI + 1, _curJ, _curK), 
            0, MPI_COMM_WORLD, &_request);
    }
}

void Process::sendDownJ(std::vector<double> downJ) {
    if (_curJ > 0) {
        MPI_Isend(downJ.data(), downJ.size(), MPI_DOUBLE, 
            getOtherRank(_curI, _curJ - 1, _curK), 
            0, MPI_COMM_WORLD, &_request);
    }
}

void Process::sendUpJ(std::vector<double> upJ) {
    if (_curJ + 1 < _countJ) {
        MPI_Isend(upJ.data(), upJ.size(), MPI_DOUBLE, 
            getOtherRank(_curI, _curJ + 1, _curK), 
            0, MPI_COMM_WORLD, &_request);
    }
}

void Process::sendDownK(std::vector<double> downK) {
    if (_curK > 0) {
        MPI_Isend(downK.data(), downK.size(), MPI_DOUBLE, 
            getOtherRank(_curI, _curJ, _curK - 1), 
            0, MPI_COMM_WORLD, &_request);
    }
}

void Process::sendUpK(std::vector<double> upK) {
    if (_curK + 1 < _countK) {
        MPI_Isend(upK.data(), upK.size(), MPI_DOUBLE, 
            getOtherRank(_curI, _curJ, _curK + 1), 
            0, MPI_COMM_WORLD, &_request);
    }
}

//Recv
std::vector<double> Process::recvDownI() {
    std::vector<double> downI (_bsize_J * _bsize_K, 0);
    if (_curI > 0) {
        MPI_Recv(downI.data(), downI.size(), MPI_DOUBLE,
            getOtherRank(_curI - 1, _curJ, _curK),
            0, MPI_COMM_WORLD, &_status);
    }
    return downI;
}

std::vector<double> Process::recvUpI() {
    std::vector<double> upI (_bsize_J * _bsize_K, 0);
    if (_curI + 1 < _countI) {
        MPI_Recv(upI.data(), upI.size(), MPI_DOUBLE,
            getOtherRank(_curI + 1, _curJ, _curK),
            0, MPI_COMM_WORLD, &_status);
    }
    return upI;
}

std::vector<double> Process::recvDownJ() {
    std::vector<double> downJ (_bsize_I * _bsize_K, 0);
    if (_curJ > 0) {
        MPI_Recv(downJ.data(), downJ.size(), MPI_DOUBLE,
            getOtherRank(_curI, _curJ - 1, _curK),
            0, MPI_COMM_WORLD, &_status);
    }
    return downJ;
}

std::vector<double> Process::recvUpJ() {
    std::vector<double> upJ (_bsize_I * _bsize_K, 0);
    if (_curJ + 1 < _countJ) {
        MPI_Recv(upJ.data(), upJ.size(), MPI_DOUBLE,
            getOtherRank(_curI, _curJ + 1, _curK),
            0, MPI_COMM_WORLD, &_status);
    }
    return upJ;
}

std::vector<double> Process::recvDownK() {
    std::vector<double> downK (_bsize_I * _bsize_J, 0);
    if (_curK > 0) {
        MPI_Recv(downK.data(), downK.size(), MPI_DOUBLE,
            getOtherRank(_curI, _curJ, _curK - 1),
            0, MPI_COMM_WORLD, &_status);
    }
    return downK;
}

std::vector<double> Process::recvUpK() {
    std::vector<double> upK (_bsize_I * _bsize_J, 0);
    if (_curK + 1 < _countK) {
        MPI_Recv(upK.data(), upK.size(), MPI_DOUBLE,
            getOtherRank(_curI, _curJ, _curK + 1),
            0, MPI_COMM_WORLD, &_status);
    }
    return upK;
}