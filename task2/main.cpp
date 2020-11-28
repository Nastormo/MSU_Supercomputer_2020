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

void init_u1(Block &b, const Block &u0, double tau, Function3D &u) {
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

void step(Block &u2, const Block& u1, const Block& u0, double tau, Function3D &u) {
    int sizeI = u2.getSizeI();
    int sizeJ = u2.getSizeJ();
    int sizeK = u2.getSizeK();

    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
                u2.getElem(i, j, k) = 2 * u1.getValElem(i, j, k) - u0.getValElem(i, j, k) + pow(tau, 2) * u1.lap_h(i, j, k);
                // if (abs(u2.getValElem(i, j, k) - u(u2.getX(i), u2.getY(j), u2.getZ(k), tau)) > 740.0) {
                //     printf("%d %d %d %d %d %d %lf\n", sizeI, sizeJ, sizeK, i, j, k, (pow(tau, 2) / 2) * u1.lap_h(i, j, k));
                //     printf("\t%lf %lf %lf\n", u1.getValElem(i - 1, j, k), u1.getValElem(i, j, k), u1.getValElem(i + 1, j, k));
                //     printf("\t%lf %lf %lf\n", u1.getValElem(i, j - 1, k), u1.getValElem(i, j, k), u1.getValElem(i, j + 1, k));
                //     printf("\t%lf %lf %lf\n", u1.getValElem(i, j, k - 1), u1.getValElem(i, j, k), u1.getValElem(i, j, k + 1));
                // }
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

    Function3D u_a(Lx, Ly, Lz);
    MPI_Request request;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    Process p(bCountI, bCountJ, bCountK, N);

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
    init_u1(massB[1], massB[0], tau, u_a);
    p.printError(massB[1], u_a, tau);
    p.update(massB[1]);

    for (int t = 2; t < K; t++) {
        step(massB[t % 3], massB[(t + 2) % 3], massB[(t + 1) % 3], tau, u_a);
        p.printError(massB[t % 3], u_a, tau * t);
        p.update(massB[t % 3]);
    }
    MPI_Finalize();
    return 0; 
}