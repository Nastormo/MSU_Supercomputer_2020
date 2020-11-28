#include "block.h"

Block::Block(int sizeI, int sizeJ, int sizeK, 
    int minI, int minJ, int minK,
    double shiftX, double shiftY, double shiftZ) 
        : _sizeI(sizeI + 2), _sizeJ(sizeJ + 2), _sizeK(sizeK + 2),
        _minI(minI), _minJ(minJ),_minK(minK),
        _shiftX(shiftX), _shiftY(shiftY), _shiftZ(shiftZ) 
{
    _raw.resize(_sizeI * _sizeJ * _sizeK, 0.0);
}

void Block::printBlock() const {
    printf("Size: %d %d %d\n", _sizeI, _sizeJ, _sizeK);
    for (int i = 1; i < _sizeI - 1; i++) {
        printf("%d:\n", i);
        for (int j = 1; j < _sizeJ - 1; j++) {
            printf("\t%d:\n", j);
            for (int k = 1; k < _sizeK - 1; k++) {
                printf("%lf ", getValElem(i, j, k));
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

void Block::printDiff(Function3D& u, double t) const {
    printf("Size: %d %d %d\n", _sizeI, _sizeJ, _sizeK);
    for (int i = 1; i < _sizeI - 1; i++) {
        printf("%d:\n", i);
        for (int j = 1; j < _sizeJ - 1; j++) {
            printf("\t%d:\n", j);
            for (int k = 1; k < _sizeK - 1; k++) {
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
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int j = 1; j < _sizeJ - 1; j++) {
            for (int k = 1; k < _sizeK - 1; k++) {
                error = std::max(abs(getValElem(i, j, k) -  u(getX(i), getY(j), getZ(k), t)), error);
            }
        }
    }
    return error;
}

double Block::getValElem(int i, int j, int k) const {
    return _raw[i * (_sizeJ * _sizeK) + j * _sizeK + k];
}

double& Block::getElem(int i, int j, int k) {
    return _raw[i * (_sizeJ * _sizeK) + j * _sizeK + k];
}

double& Block::operator()(int i, int j, int k) {
    return _raw[(i + 1 - _minI) * (_sizeJ * _sizeK) + 
        (j + 1 - _minJ) * _sizeK + (k + 1 - _minK)];
}  

//TODO
double Block::operator()(int i, int j, int k) const {
    return _raw[i * (_sizeJ * _sizeK) + j * _sizeK + k];
}

//Uncorrect 0 -> 1 size* - 1 -> size* - 1
//GetX
std::vector<double> Block::getDownI() const { 
    std::vector<double> downX ((_sizeJ - 2) * (_sizeK - 2), 0);
    for (int j = 1; j < _sizeJ - 1; j++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            downX[(j - 1) * (_sizeK - 2) + (k - 1)] = getValElem(1, j, k);
        }
    }
    return downX;
}

std::vector<double> Block::getUpI() const { 
    std::vector<double> upX ((_sizeJ - 2) * (_sizeK - 2), 0);
    for (int j = 1; j < _sizeJ - 1; j++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            upX[(j - 1) * (_sizeK - 2) + (k - 1)] = getValElem(_sizeI - 2, j, k);
        }
    }
    return upX;
}

//GetY
std::vector<double> Block::getDownJ() const { 
    std::vector<double> downY ((_sizeI - 2) * (_sizeK - 2), 0);
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            downY[(i - 1) *  (_sizeK - 2) + (k - 1)] = getValElem(i, 1, k);
        }
    }
    return downY;
}

std::vector<double> Block::getUpJ() const { 
    std::vector<double> upY ((_sizeI - 2) * (_sizeK - 2), 0);
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            upY[(i - 1) * (_sizeK - 2) + (k - 1)] = getValElem(i, _sizeJ - 2, k);
        }
    }
    return upY;
}

//GetZ
std::vector<double> Block::getDownK() const { 
    std::vector<double> downZ ((_sizeI - 2) * (_sizeJ - 2), 0);
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int j = 1; j < _sizeJ - 1; j++) {
            downZ[(i - 1) * (_sizeJ - 2) + (j - 1)] = getValElem(i, j, 1);
        }
    }
    return downZ;
}

std::vector<double> Block::getUpK() const { 
    std::vector<double> upZ ((_sizeI - 2) * (_sizeJ - 2), 0);
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int j = 1; j < _sizeJ - 1; j++) {
            upZ[(i - 1) * (_sizeJ - 2) + (j - 1)] = getValElem(i, j, _sizeK - 2);
        }
    }
    return upZ;
}

//SetX
void Block::setDownI(const std::vector<double>& downI) {
    for (int j = 1; j < _sizeJ - 1; j++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            getElem(0, j, k) = downI[(j - 1) * (_sizeK - 2) + (k - 1)];
        }
    }
}

void Block::setUpI(const std::vector<double>& upI) {
    for (int j = 1; j < _sizeJ - 1; j++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            getElem(_sizeI - 1, j, k) = upI[(j - 1) * (_sizeK - 2) + (k - 1)];
        }
    }
}

//SetY
void Block::setDownJ(const std::vector<double>& downJ) {
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            getElem(i, 0, k) = downJ[(i - 1) * (_sizeK - 2) + (k - 1)];
        }
    }
}

void Block::setUpJ(const std::vector<double>& upJ) {
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            getElem(i, _sizeJ - 1, k) = upJ[(i - 1) * (_sizeK - 2) + (k - 1)];
        }
    }
}

//SetZ
void Block::setDownK(const std::vector<double>& downK) {
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int j = 1; j < _sizeJ - 1; j++) {
            getElem(i, j, 0) = downK[(i - 1) * (_sizeJ - 2) + (j - 1)];
        }
    }
}

void Block::setUpK(const std::vector<double>& upK) {
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int j = 1; j < _sizeJ - 1; j++) {
            getElem(i, j, _sizeK - 1) = upK[(i - 1) * (_sizeJ - 2) + (j - 1)];
        }
    }
}

double Block::lap_h(int i, int j, int k) const {
    return (getValElem(i - 1, j, k) - 2 * getValElem(i, j, k) + getValElem(i + 1, j, k)) / pow(_shiftX, 2) + 
        (getValElem(i, j - 1, k) - 2 * getValElem(i, j, k) + getValElem(i, j + 1, k)) / pow(_shiftY, 2) +
        (getValElem(i, j, k - 1) - 2 * getValElem(i, j, k) + getValElem(i, j, k + 1)) / pow(_shiftZ, 2); 
}
