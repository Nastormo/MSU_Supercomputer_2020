#include "parallel.h"


Parallel::Parallel(Config& conf, Function3D& u)
    : _conf(conf), _p(conf), _u(u), _t(0), _tau(conf.getTau()), _K(conf.getK())
{
    std::vector<double>L = _u.getL();
    int N = _conf.getN();

    std::vector<int> bsize = _p.getSize();
    std::vector<int> bmin = _p.getMin();
    std::vector<double> bshift = { L[0]/(N + 1), L[1]/(N + 1), L[2]/(N + 1) };
    
    for (int i = 0; i < 3; i++) {
        _massB.push_back(Block(bsize, bmin, bshift));
    }
    _bshift = _massB[0].getShift();
    _bsize = _massB[0].getSize();
    _bmin = _massB[0].getMin();
}

void Parallel::init_u0() {
    Block& u0 = _massB[_t % 3];

    #pragma omp parallel for
    for (int i = 1; i < _bsize[0] - 1; i++) {
        for (int j = 1; j < _bsize[1] - 1; j++) {
            for (int k = 1; k < _bsize[2] - 1; k++) {
                u0.getElem(i, j, k) = _u(u0.getX(i), u0.getY(j), u0.getZ(k), 0);
            }
        }
    }
}

void Parallel::init_u1() {
    Block& u0 = _massB[(_t - 1) % 3];
    Block& u1 = _massB[_t % 3];

    #pragma omp parallel for
    for (int i = 1; i < _bsize[0] - 1; i++) {
        for (int j = 1; j < _bsize[1] - 1; j++) {
            for (int k = 1; k < _bsize[2] - 1; k++) {
                u1.getElem(i, j, k) = u0.getValElem(i, j, k) + (pow(_tau, 2) / 2) * u0.lap_h(i, j, k);
            }
        }
    } 
}

void Parallel::step() {
    Block& u0 = _massB[(_t - 2) % 3];
    Block& u1 = _massB[(_t - 1) % 3];
    Block& u2 = _massB[_t % 3];

    #pragma omp parallel for
    for (int i = 1; i < _bsize[0] - 1; i++) {
        for (int j = 1; j < _bsize[1] - 1; j++) {
            for (int k = 1; k < _bsize[2] - 1; k++) {
                u2.getElem(i, j, k) = 2 * u1.getValElem(i, j, k) - u0.getValElem(i, j, k) + pow(_tau, 2) * u1.lap_h(i, j, k);
            }
        }
    }
}

void Parallel::printError() {
    double time = _t * _tau;
    _p.printError(time, _massB[_t % 3].getError(_u, time));
}

void Parallel::update() {
    _p.update(_massB[_t % 3]);
}

void Parallel::process() {
    init_u0();
    printError();
    update();
    _t++;
    
    init_u1();
    printError();
    update();
    _t++;

    while (_t < _K) {
        step();
        printError();
        update();
        _t++;
    }
}