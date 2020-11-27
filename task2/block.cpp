#include "block.h"

Block::Block(int sizeI, int sizeJ, int sizeK, 
    int minI, int minJ, int minK,
    double shiftX, double shiftY, double shiftZ) 
        : _sizeI(sizeI + 2), _sizeJ(sizeJ + 2), _sizeK(sizeK + 2),
        _minI(minI), _minJ(minJ),_minK(minK),
        _shiftX(shiftX), _shiftY(shiftY), _shiftZ(shiftZ) 
{
    _raw.resize(sizeI * sizeJ * sizeK, 0.0);
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
std::vector<double> Block::getDownI() { 
    std::vector<double> downX ((_sizeJ - 2) * (_sizeK - 2), 0);
    for (int j = 1; j < _sizeJ - 1; j++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            downX[j * (_sizeK - 2) + k] = this->getValElem(1, j, k);
        }
    }
    return downX;
}

std::vector<double> Block::getUpI() { 
    std::vector<double> upX ((_sizeJ - 2) * (_sizeK - 2), 0);
    for (int j = 1; j < _sizeJ - 1; j++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            upX[j * (_sizeK - 2) + k] = this->getValElem(_sizeI - 2, j, k);
        }
    }
    return upX;
}

//GetY
std::vector<double> Block::getDownJ() { 
    std::vector<double> downY ((_sizeI - 2) * (_sizeK - 2), 0);
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            downY[i * (_sizeJ * _sizeK) + k] = this->getValElem(i, 1, k);
        }
    }
    return downY;
}

std::vector<double> Block::getUpJ() { 
    std::vector<double> upY ((_sizeI - 2) * (_sizeK - 2), 0);
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            upY[i * (_sizeJ * _sizeK) + k] = this->getValElem(i, _sizeJ - 2, k);
        }
    }
    return upY;
}

//GetZ
std::vector<double> Block::getDownK() { 
    std::vector<double> downZ ((_sizeI - 2) * (_sizeJ - 2), 0);
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int j = 1; j < _sizeJ - 1; j++) {
            downZ[i * (_sizeJ * _sizeK) + j * _sizeK] = this->getValElem(i, j, 1);
        }
    }
    return downZ;
}

std::vector<double> Block::getUpK() { 
    std::vector<double> upZ ((_sizeI - 2) * (_sizeJ - 2), 0);
    for (int i = 1; i < _sizeI - 1; i++) {
        for (int j = 1; j < _sizeJ - 1; j++) {
            upZ[i * (_sizeJ * _sizeK) + j * _sizeK] = this->getValElem(i, j, _sizeK - 2);
        }
    }
    return upZ;
}

//SetX
void Block::setDownX(std::vector<double>& downX) {
    for (int j = 0; j < _sizeJ; j++) {
        for (int k = 0; k < _sizeK; k++) {
            this->getElem(0, j, k) = downX[j * _sizeK + k];
        }
    }
}

void Block::setUpX(std::vector<double>& upX) {
    for (int j = 0; j < _sizeJ; j++) {
        for (int k = 0; k < _sizeK; k++) {
            this->getElem(_sizeI - 1, j, k) = upX[j * _sizeK + k];
        }
    }
}

//SetY
void Block::setDownY(std::vector<double>& downY) {
    for (int i = 0; i < _sizeI; i++) {
        for (int k = 0; k < _sizeK; k++) {
            this->getElem(i, 0, k) = downY[i * (_sizeJ * _sizeK) + k];
        }
    }
}

void Block::setUpY(std::vector<double>& upY) {
    for (int i = 0; i < _sizeI; i++) {
        for (int k = 0; k < _sizeK; k++) {
            this->getElem(i, _sizeJ - 1, k) = upY[i * (_sizeJ * _sizeK) + k];
        }
    }
}

//SetZ
void Block::setDownZ(std::vector<double>& downZ) {
    for (int i = 0; i < _sizeI; i++) {
        for (int j = 0; j < _sizeJ; j++) {
            this->getElem(i, j, 0) = downZ[i * (_sizeJ * _sizeK) + j * _sizeK];
        }
    }
}

void Block::setUpZ(std::vector<double>& upZ) {
    for (int i = 0; i < _sizeI; i++) {
        for (int j = 0; j < _sizeJ; j++) {
            this->getElem(i, j, _sizeK - 1) = upZ[i * (_sizeJ * _sizeK) + j * _sizeK];
        }
    }
}

double Block::lap_h(int i, int j, int k) const {
    return (getValElem(i - 1, j, k) - 2 * getValElem(i, j, k) + getValElem(i + 1, j, k)) / pow(_shiftX, 2) + 
        (getValElem(i, j - 1, k) - 2 * getValElem(i, j, k) + getValElem(i, j + 1, k)) / pow(_shiftY, 2) +
        (getValElem(i, j, k - 1) - 2 * getValElem(i, j, k) + getValElem(i, j, k + 1)) / pow(_shiftZ, 2); 
}
