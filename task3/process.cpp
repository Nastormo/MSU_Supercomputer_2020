#include "process.h"
#include "parallel.hpp"

Process::Process(Config conf)
{
    int N = conf.getN();
    int argc = conf.getArgc();
    char** argv = conf.getArgv();

    _ndims.resize(3, 0);
    _position.resize(3, 0);
    _bsize.resize(3, 0);
    _bmin.resize(3, 0);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &_p_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

    MPI_Dims_create(_p_count, 3, _ndims.data());

    conf.setNdims(_ndims.data());
    if (_rank == 0) {
        conf.printConfig();
    }

    _position[0] = _rank / (_ndims[1] * _ndims[2]);
    _position[1] = (_rank - (_ndims[1] * _ndims[2]) * _position[0]) / _ndims[2];
    _position[2] = _rank % _ndims[2];

    for (int i = 0; i < 3; i++) {
        int blockSize = (N + _ndims[i] - 1) / _ndims[i];
        _bsize[i] = std::min(N, 
            blockSize * (_position[i] + 1)) - blockSize * _position[i];
        _bmin[i] = blockSize * _position[i];
    }

    _startTime = MPI_Wtime();
}

Process::~Process() {
    if(_rank == 0) {
        printf("Elapsed time = %lf\n", MPI_Wtime() - _startTime);
    }
    MPI_Finalize();
}

void Process::printError(const Block& b, Function3D &u, double t) {
    double otherError = getError(b, u, t);
    double error; 

    MPI_Reduce(&otherError, &error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (_rank == 0) {
        std::cout << "t: " << t << " Error: " << error << std::endl;
    }
}

void Process::update(Block &b) {
    std::vector<double> rUSlice;
    std::vector<double> rDSlice;
    std::vector<double> sUSlice;
    std::vector<double> sDSlice;
    std::vector<int> size = b.getSize();

    for (int axis = 0; axis < 3; axis++) {
        sDSlice = b.getSlice(axis, 1);
        sUSlice = b.getSlice(axis, size[axis] - 2);

        //Send down, recv up
        Send(sDSlice, axis, -1);
        rUSlice = Recv(axis, 1);
        b.setSlice(rUSlice, axis, size[axis] - 1);

        //Send up, recv down
        Send(sUSlice, axis, 1);
        rDSlice = Recv(axis, -1);
        b.setSlice(rDSlice, axis, 0);

        // if(_rank == 1) std::cout << axis << ' ' << rDSlice[10] << std::endl;
        // if(_rank == 0) std::cout << axis << ' ' << sUSlice[10] << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

//shift only 1 or -1
void Process::Send(std::vector<double>& slice, int axis, int shift) {
    if (_position[axis] + shift >= 0 && 
        _position[axis] + shift < _ndims[axis]) {
        std::vector<int> new_position(_position);
        new_position[axis] += shift;
        MPI_Isend(slice.data(), slice.size(), MPI_DOUBLE, 
            getOtherRank(new_position), 
            0, MPI_COMM_WORLD, &_request);
    }
}

//shift only 1 or -1
std::vector<double> Process::Recv(int axis, int shift) {
    std::vector<double> slice (
        _bsize[(axis + 1) % 3] * _bsize[(axis + 2) % 3], 0
    );
    if(_position[axis] + shift >= 0 &&
        _position[axis] + shift < _ndims[axis]) {
        std::vector<int> new_position(_position);
        new_position[axis] += shift;
        MPI_Recv(slice.data(), slice.size(), MPI_DOUBLE,
            getOtherRank(new_position),
            0, MPI_COMM_WORLD, &_status);
    }

    return slice;
}

int Process::getOtherRank(const std::vector<int> &position) {
    return position[0] * (_ndims[1] * _ndims[2]) + 
        position[1] * _ndims[2] + 
        position[2]; 
}