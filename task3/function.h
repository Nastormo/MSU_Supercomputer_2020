#pragma once
#include <math.h>
#include <stdio.h>
#include <algorithm>

class Function3D {
    double _Lx;
    double _Ly;
    double _Lz;
public:
    Function3D(double Lx, double Ly, double Lz);
    double getLx() { return _Lx; }
    double getLy() { return _Ly; }
    double getLz() { return _Lz; }
    double a_t();
    double operator()(double x, double y, double z, double t);
};