#include "block.h"

Block::Block(int sizeI, int sizeJ, int sizeK, int minI, int minJ, int minK) 
        : _sizeI(sizeI + 2), _sizeJ(sizeJ + 2), _sizeK(sizeK + 2),
        _minI(minI), _minJ(minJ),_minK(minK) 
{
    _raw.resize(sizeI * sizeJ * sizeK, 0.0);
}

double Block::getValElem(int i, int j, int k) {
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

//GetX
std::vector<double> Block::getDownX() { 
    std::vector<double> downX ((_sizeJ - 2) * (_sizeK - 2), 0);
    for (int j = 1; j < _sizeJ - 1; j++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            downX[j * (_sizeK - 2) + k] = this->getValElem(0, j, k);
        }
    }
    return downX;
}

std::vector<double> Block::getUpX() { 
    std::vector<double> upX ((_sizeJ - 2) * (_sizeK - 2), 0);
    for (int j = 1; j < _sizeJ - 1; j++) {
        for (int k = 1; k < _sizeK - 1; k++) {
            upX[j * (_sizeK - 2) + k] = this->getValElem(_sizeI - 1, j, k);
        }
    }
    return upX;
}

//GetY
std::vector<double> Block::getDownY() { 
    std::vector<double> downY (_sizeI * _sizeK, 0);
    for (int i = 0; i < _sizeI; i++) {
        for (int k = 0; k < _sizeK; k++) {
            downY[i * (_sizeJ * _sizeK) + k] = this->getValElem(i, 0, k);
        }
    }
    return downY;
}

std::vector<double> Block::getUpY() { 
    std::vector<double> upY (_sizeI * _sizeK, 0);
    for (int i = 0; i < _sizeI; i++) {
        for (int k = 0; k < _sizeK; k++) {
            upY[i * (_sizeJ * _sizeK) + k] = this->getValElem(i, _sizeJ - 1, k);
        }
    }
    return upY;
}

//GetZ
std::vector<double> Block::getDownZ() { 
    std::vector<double> downZ (_sizeI * _sizeJ, 0);
    for (int i = 0; i < _sizeI; i++) {
        for (int j = 0; j < _sizeJ; j++) {
            downZ[i * (_sizeJ * _sizeK) + j * _sizeK] = this->getValElem(i, j, 0);
        }
    }
    return downZ;
}

std::vector<double> Block::getUpZ() { 
    std::vector<double> upZ (_sizeI * _sizeJ, 0);
    for (int i = 0; i < _sizeI; i++) {
        for (int j = 0; j < _sizeJ; j++) {
            upZ[i * (_sizeJ * _sizeK) + j * _sizeK] = this->getValElem(i, j, _sizeK - 1);
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
