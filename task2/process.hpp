#include <mpi.h>

class Process {
    int _rank;
    int _p_count;

    int _countX;
    int _countY;
    int _countZ;

    int _curX;
    int _curY;
    int _curZ;

    int _bsize_X;
    int _bsize_Y;
    int _bsize_Z;

    int _bmin_i;
    int _bmin_j;
    int _bmin_k;

public:  
    Process(int countX, int countY, int countZ, int N) 
            : _countX(countX),
            _countY(countY),
            _countZ(countZ) {
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &_p_count);

        _curX = _rank / (_countY * _countZ);;
        _curY = (_rank - (_countY * _countZ) * _curX) / _countZ;
        _curZ = _rank % _countZ;

        _bsize_X = N / _countX;
        _bsize_Y = N / _countY;
        _bsize_Z = N / _countZ;

        _bmin_i = _bsize_X * _curX;
        _bmin_j = _bsize_Y * _curY;
        _bmin_k = _bsize_Z * _curZ;
    }

    int getRank() { return _rank; }

    int getX() { return _curX; }
    int getY() { return _curY; }
    int getZ() { return _curZ; }

    int getSizeX() { return _bsize_X; }
    int getSizeY() { return _bsize_Y; }
    int getSizeZ() { return _bsize_Z; }

    int get_i() { return _bmin_i; }
    int get_j() { return _bmin_j; }
    int get_k() { return _bmin_k; }
};