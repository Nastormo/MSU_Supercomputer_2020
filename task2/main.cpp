#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#include "function.h"
#include "block.h"
#include "process.hpp"


void init_u0(Block &b, Function3D &u) {
    int startI = b.getMinI() - 1;
    int startJ = b.getMinJ() - 1;
    int startK = b.getMinK() - 1;

    int sizeX = b.getSizeI();
    int sizeY = b.getSizeJ();
    int sizeZ = b.getSizeK();

    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            for (int z = 0; z < sizeZ; z++) {
                b.getElem(x, y, z) = u(0, 0, 0, 0);
            }
        }
    }

}

int main(int argc, char** argv) {
    if (argc < 9) {
        std::cout << "A few argument" << std::endl;
        return 1;
    }
    double L_x = atof(argv[1]);
    double L_y = atof(argv[2]);
    double L_z = atof(argv[3]);
    int N = atoi(argv[4]);
    int K = atoi(argv[5]);
    int bCountX = atoi(argv[6]);
    int bCountY = atoi(argv[7]);
    int bCountZ = atoi(argv[8]);

    Function3D u_a(L_x, L_y, L_z);

    MPI_Init(&argc, &argv);


    Process p(bCountX, bCountY, bCountZ, N);

    std::vector<Block> massB (3, 
        Block(p.getSizeX(), p.getSizeY(), p.getSizeZ(), p.get_i(), p.get_j(), p.get_k())
    );
    init_u0(massB[0], u_a);
    //std::cout << p.getRank();
    
    MPI_Finalize();
    return 0; 
}