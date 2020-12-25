#pragma once

#include "config.h"
#include "process.h"
#include "function.h"
#include "block.h"

class Parallel {
    Config _conf;
    Process _p;
    Function3D _u;

    std::vector<double> _bshift;
    std::vector<int> _bsize;
    std::vector<int> _bmin;

    int _t;
    int _K;
    double _tau;

    std::vector<Block> _massB; 
public:
    Parallel(Config& conf, Function3D& u);
    void init_u0();
    void init_u1();
    void step();
    void printError();
    void update();
    void process();
};