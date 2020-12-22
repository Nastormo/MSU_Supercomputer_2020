#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <sys/stat.h> 
#include <sys/types.h> 


// #include "block.h"

#include "config.h"
#include "function.h"
#include "process.h"
#include "omp.h"
#include "cuda.hpp"

// void saveMatrix(Block &b, Function3D &u, Process &p, double t, int N) {
//     int sizeI = b.getSizeI();
//     int sizeJ = b.getSizeJ();
//     int sizeK = b.getSizeK();
//     std::vector<double> vb((sizeI - 2) * (sizeJ - 2) * (sizeK - 2));
//     std::vector<double> vu((sizeI - 2) * (sizeJ - 2) * (sizeK - 2));

//     for (int i = 1; i < sizeI - 1; i++) {
//         for (int j = 1; j < sizeJ - 1; j++) {
//             for (int k = 1; k < sizeK - 1; k++) {
//                 vb[(i - 1) * ((sizeJ - 2) * (sizeK - 2)) + (j - 1) * (sizeK - 2) + (k - 1)] = b.getValElem(i, j, k);
//                 vu[(i - 1) * ((sizeJ - 2) * (sizeK - 2)) + (j - 1) * (sizeK - 2) + (k - 1)] = u(b.getX(i), b.getY(j), b.getZ(k), t);
//             }
//         }
//     }
//     std::stringstream path;
//     path << "matrix/" << N;
//     mkdir(path.str().c_str(), 0777);
//     path << "/" << p.getPCount();
//     mkdir(path.str().c_str(), 0777);
//     path << "/" << p.getI() << p.getJ() << p.getK();
//     mkdir(path.str().c_str(), 0777);
//     std::stringstream name;
//     name << "_" << t;
//     std::string nb = path.str() + "/B" + name.str();
//     std::ofstream outFileB(nb.c_str(), std::ios::out | std::ios::binary);
//     outFileB.write((char *)vb.data(), vb.size() * sizeof(double));

//     std::string nu = path.str() + "/U" + name.str();
//     std::ofstream outFileU(nu.c_str(), std::ios::out | std::ios::binary);
//     outFileU.write((char *)vu.data(), vu.size() * sizeof(double));
// }



int main(int argc, char** argv) {
    Config conf(argc, argv);

    test();
    return 0;

    // Function3D u_a(conf.getLx(), conf.getLy(), conf.getLz());
    // Process p(conf);

    // double Lx = conf.getLx();
    // double Ly = conf.getLy();
    // double Lz = conf.getLz();
    // double tau = conf.getTau();
    // int K = conf.getK();
    // int N = conf.getN();
    // std::vector<int> bsize = p.getSize();
    // std::vector<int> bmin = p.getMin();
    // std::vector<double> bshift = { Lx/(N + 1), Ly/(N + 1), Lz/(N + 1) };

    // std::vector<Block> massB (3,
    //     Block(bsize, bmin, bshift)
    // );

    // init_u0(massB[0], u_a);
    // p.printError(massB[0], u_a, 0.0);
    // p.update(massB[0]);
    // //saveMatrix(massB[0], u_a, p, 0, N);
    
    // init_u1(massB[1], massB[0], tau, u_a);
    // p.printError(massB[1], u_a, tau);
    // p.update(massB[1]);
    // // saveMatrix(massB[1], u_a, p, tau, N);

    // std::cout << K << std::endl;
    // for (int t = 2; t < K; t++) {
    //     step(massB[t % 3], massB[(t + 2) % 3], massB[(t + 1) % 3], tau, u_a);
    //     p.printError(massB[t % 3], u_a, tau * t);
    //     p.update(massB[t % 3]);
    //     // saveMatrix(massB[t % 3], u_a, p, tau * t, N);
    // }
    
    // return 0; 
}