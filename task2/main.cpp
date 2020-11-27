#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#include "function.h"
#include "block.h"
#include "process.h"


void init_u0(Block &b, Function3D &u) {
    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();

    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
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
    int bCountI = atoi(argv[7]);
    int bCountJ = atoi(argv[8]);
    int bCountK = atoi(argv[9]);
    double tau = T / K;

    //printf("Lx: %lf; Ly: %lf; Lz: %lf;\nT: %lf; N: %d; K: %d;\nCountI: %d; CountJ: %d; CountK: %d;\n", Lx, Ly, Lz, T, N, K, bCountI, bCountJ, bCountK);

    Function3D u_a(Lx, Ly, Lz);
    MPI_Request request;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    Process p(bCountI, bCountJ, bCountK, N, request, status);

    // std::vector<Block> massB;
    // for (int i = 0; i < 3; i++) {
    //     massB.push_back( 
    //         Block(
    //             p.getSizeX(), p.getSizeY(), p.getSizeZ(),
    //             p.get_i(), p.get_j(), p.get_k(),
    //             Lx / N, Ly / N, Lz / N
    //         )
    //     );
    // }
    std::vector<Block> massB (3,
        Block(
            p.getSizeX(), p.getSizeY(), p.getSizeZ(),
            p.get_i(), p.get_j(), p.get_k(),
            Lx / N, Ly / N, Lz / N
        )
    );

    init_u0(massB[0], u_a);
    p.printError(massB[0], u_a, 0.0);
    p.update(massB[0]);
    init_u1(massB[1], massB[0], tau);
    p.printError(massB[1], u_a, tau);
    p.update(massB[1]);

    //p.sendDownI(massB[1].getDownI(), request);
    //p.sendUpI(massB[1].getUpI(), request);
    //p.sendDownJ(massB[1].getDownJ(), request);
    //p.sendUpJ(massB[1].getUpJ(), request);
    //p.sendDownK(massB[1].getDownK(), request);
    //p.sendUpK(massB[1].getUpK(), request);

    // for (int t = 1; t < K; t++) {


    // }
    //std::cout << p.getRank();
    //massB.clear();
    MPI_Finalize();
    return 0; 
}