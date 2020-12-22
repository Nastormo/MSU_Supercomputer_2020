#pragma once

#include "function.h"
#include "block.h"

void test();
void init_u0(Block &b, Function3D &u);
void init_u1(Block &b, const Block &u0, double tau, Function3D &u);
void step(Block &u2, const Block& u1, const Block& u0, double tau, Function3D &u);