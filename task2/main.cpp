#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#include "function.h"
#include "block.h"
#include "process.hpp"


void init_u0(Block &b, Function3D &u) {
    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();

    for (int i = 0; i < sizeI; i++) {
        for (int j = 0; j < sizeJ; j++) {
            for (int k = 0; k < sizeK; k++) {
                b.getElem(i, j, k) = u(b.getX(i), b.getY(j), b.getZ(k), 0);
            }
        }
    }
}

void init_u1(Block &b, const Block &u0, double tau) {
    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();

    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
                b.getElem(i, j, k) = u0.getValElem(i, j, k) + (pow(tau, 2) / 2) * u0.lap_h(i, j, k);
            }
        }
    } 
}

int main(int argc, char** argv) {
    if (argc < 9) {
        std::cout << "A few argument" << std::endl;
        return 1;
    }
    double Lx = atof(argv[1]);
    double Ly = atof(argv[2]);
    double Lz = atof(argv[3]);
    double T = atof(argv[4]);
    int N = atoi(argv[5]);
    int K = atoi(argv[6]);
    int bCountX = atoi(argv[7]);
    int bCountY = atoi(argv[8]);
    int bCountZ = atoi(argv[9]);

    Function3D u_a(Lx, Ly, Lz);

    MPI_Init(&argc, &argv);

    Process p(bCountX, bCountY, bCountZ, N);

    std::vector<Block> massB (3, 
        Block(
            p.getSizeX(), p.getSizeY(), p.getSizeZ(),
            p.get_i(), p.get_j(), p.get_k(),
            Lx / N, Ly / N, Lz / N
        )
    );

    init_u0(massB[0], u_a);
    init_u1(massB[1], massB[0], T / K);
    //std::cout << p.getRank();
    
    MPI_Finalize();
    return 0; 
}