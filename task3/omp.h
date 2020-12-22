#pragma once

#include "block.h"
#include "function.h"

class Omp {

public:
    void init_u0(Block &b, Function3D &u);
    void init_u1(Block &b, const Block &u0, double tau, Function3D &u);
    void step(Block &u2, const Block& u1, const Block& u0, double tau, Function3D &u);
};