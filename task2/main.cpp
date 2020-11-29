#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <sys/stat.h> 
#include <sys/types.h> 

#include "function.h"
#include "block.h"
#include "process.h"

void saveMatrix(Block &b, Function3D &u, Process &p, double t, int N) {
    int sizeI = b.getSizeI();
    int sizeJ = b.getSizeJ();
    int sizeK = b.getSizeK();
    std::vector<double> vb((sizeI - 2) * (sizeJ - 2) * (sizeK - 2));
    std::vector<double> vu((sizeI - 2) * (sizeJ - 2) * (sizeK - 2));

    for (int i = 1; i < sizeI - 1; i++) {
        for (int j = 1; j < sizeJ - 1; j++) {
            for (int k = 1; k < sizeK - 1; k++) {
                vb[(i - 1) * ((sizeJ - 2) * (sizeK - 2)) + (j - 1) * (sizeK - 2) + (k - 1)] = b.getValElem(i, j, k);
                vu[(i - 1) * ((sizeJ - 2) * (sizeK - 2)) + (j - 1) * (sizeK - 2) + (k - 1)] = u(b.getX(i), b.getY(j), b.getZ(k), t);
            }
        }
    }
    std::stringstream path;
    path << "matrix/" << N;
    mkdir(path.str().c_str(), 0777);
    path << "/" << p.getI() << p.getJ() << p.getK();
    mkdir(path.str().c_str(), 0777);
    std::stringstream name;
    name << "_" << t;
    std::string nb = path.str() + "/B" + name.str();
    std::ofstream outFileB(nb.c_str(), std::ios::out | std::ios::binary);
    outFileB.write((char *)vb.data(), vb.size() * sizeof(double));

    std::string nu = path.str() + "/U" + name.str();
    std::ofstream outFileU(nu.c_str(), std::ios::out | std::ios::binary);
    outFileU.write((char *)vu.data(), vu.size() * sizeof(double));
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

    Process p(bCountI, bCountJ, bCountK, N, argc, argv);

    std::vector<Block> massB (3,
        Block(
            p.getSizeX(), p.getSizeY(), p.getSizeZ(),
            p.get_i(), p.get_j(), p.get_k(),
            Lx / (N + 1), Ly / (N + 1), Lz / (N + 1)
        )
    );

    init_u0(massB[0], u_a);
    p.printError(massB[0], u_a, 0.0);
    p.update(massB[0]);
    //saveMatrix(massB[0], u_a, p, 0, N);

    init_u1(massB[1], massB[0], tau, u_a);
    p.printError(massB[1], u_a, tau);
    p.update(massB[1]);
    //saveMatrix(massB[1], u_a, p, tau, N);

    for (int t = 2; t < K; t++) {
        step(massB[t % 3], massB[(t + 2) % 3], massB[(t + 1) % 3], tau, u_a);
        p.printError(massB[t % 3], u_a, tau * t);
        p.update(massB[t % 3]);
        //saveMatrix(massB[t % 3], u_a, p, tau * t, N);
    }
    
    return 0; 
}