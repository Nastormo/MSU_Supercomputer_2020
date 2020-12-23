#include "block.h"
#include <iostream>

Block::Block(std::vector<int>& size, std::vector<int>& min,
    std::vector<double>& shift) 
        : _size(size), _min(min), _shift(shift)
{
    for (int i = 0; i < _size.size(); i++) {
        _size[i] += 2;
    }
    _raw.resize(_size[0] * _size[1] * _size[2], 0.0);
}

void Block::printBlock() const {
    printf("Size: %d %d %d\n", _size[0], _size[1], _size[2]);
    for (int i = 1; i < _size[0] - 1; i++) {
        printf("%d:\n", i);
        for (int j = 1; j < _size[1] - 1; j++) {
            printf("\t%d:\n", j);
            for (int k = 1; k < _size[2] - 1; k++) {
                printf("%lf ", getValElem(i, j, k));
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

void Block::saveBlock(std::string &str) const {
    std::ofstream outFile(str.c_str(), std::ios::out | std::ios::binary);
    outFile.write((char *)_raw.data(), _raw.size() * sizeof(double));
}

void Block::printDiff(Function3D& u, double t) const {
    printf("Size: %d %d %d\n", _size[0], _size[1], _size[2]);
    for (int i = 1; i < _size[0] - 1; i++) {
        printf("%d:\n", i);
        for (int j = 1; j < _size[1] - 1; j++) {
            printf("\t%d:\n", j);
            for (int k = 1; k < _size[2] - 1; k++) {
                printf("%lf ", getValElem(i, j, k) - u(getX(i), getY(j), getZ(k), t));
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

double Block::getError(Function3D& u, double t) const {
    double error = 0;
    #pragma omp parallel for
    for (int i = 1; i < _size[0] - 1; i++) {
        for (int j = 1; j < _size[1] - 1; j++) {
            for (int k = 1; k < _size[2] - 1; k++) {
                // #pragma omp critical 
                // { 
                error = std::max(std::abs(getValElem(i, j, k) -  u(getX(i), getY(j), getZ(k), t)), error);
                //}
            }
        }
    }
    return error;
}

double Block::getValElem(int i, int j, int k) const {
    return _raw[i * (_size[1] * _size[2]) + j * _size[2] + k];
}

double& Block::getElem(int i, int j, int k) {
    return _raw[i * (_size[1] * _size[2]) + j * _size[2] + k];
}

double& Block::operator()(int i, int j, int k) {
    return _raw[(i + 1 - _min[0]) * (_size[1] * _size[2]) + 
        (j + 1 - _min[1]) * _size[2] + (k + 1 - _min[2])];
}  

//TODO
double Block::operator()(int i, int j, int k) const {
    return _raw[i * (_size[1] * _size[2]) + j * _size[2] + k];
}

//GetX
std::vector<double> Block::getSlice(int axis, int item) const {
    int i_min[3] = {1, 1, 1};
	int i_max[3] = {_size[0] - 1, _size[1] - 1, _size[2] - 1};

    i_min[axis] = item;
    i_max[axis] = item + 1;

    std::vector<double> slice;

    for (int i = i_min[0]; i < i_max[0]; i++) {
        for (int j = i_min[1]; j < i_max[1]; j++) {
            for (int k = i_min[2]; k < i_max[2]; k++) {
                slice.push_back(getValElem(i, j, k));
            }
        }
    }
    return slice;
}

//SetX
void Block::setSlice(const std::vector<double>& slice, int axis, int item) {
    int i_min[3] = {1, 1, 1};
	int i_max[3] = {_size[0] - 1, _size[1] - 1, _size[2] - 1};

    i_min[axis] = item;
    i_max[axis] = item + 1;

    int ind = 0;
    for (int i = i_min[0]; i < i_max[0]; i++) {
        for (int j = i_min[1]; j < i_max[1]; j++) {
            for (int k = i_min[2]; k < i_max[2]; k++) {
                getElem(i, j, k) = slice.at(ind++);
            }
        }
    }
}

double Block::lap_h(int i, int j, int k) const {
    return (getValElem(i - 1, j, k) - 2 * getValElem(i, j, k) + getValElem(i + 1, j, k)) / pow(_shift[0], 2) + 
        (getValElem(i, j - 1, k) - 2 * getValElem(i, j, k) + getValElem(i, j + 1, k)) / pow(_shift[1], 2) +
        (getValElem(i, j, k - 1) - 2 * getValElem(i, j, k) + getValElem(i, j, k + 1)) / pow(_shift[2], 2); 
}
