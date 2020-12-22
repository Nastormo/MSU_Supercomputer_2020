#include "parallel.hpp"

void test() {
    return;
}

void init_u0(Block &b, Function3D &u) {
    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();

    #pragma omp parallel for
    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
                b.getElem(i, j, k) = u(b.getX(i), b.getY(j), b.getZ(k), 0);
            }
        }
    }
}

void init_u1(Block &b, const Block &u0, double tau, Function3D &u) {
    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();

    #pragma omp parallel for
    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
                b.getElem(i, j, k) = u0.getValElem(i, j, k) + (pow(tau, 2) / 2) * u0.lap_h(i, j, k);
            }
        }
    } 
}

void step(Block &u2, const Block& u1, const Block& u0, double tau, Function3D &u) {
    int sizeI = u2.getSizeI();
    int sizeJ = u2.getSizeJ();
    int sizeK = u2.getSizeK();

    #pragma omp parallel for
    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
                u2.getElem(i, j, k) = 2 * u1.getValElem(i, j, k) - u0.getValElem(i, j, k) + pow(tau, 2) * u1.lap_h(i, j, k);
            }
        }
    }
}