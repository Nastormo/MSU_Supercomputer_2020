#pragma once
#include <math.h>

class Function3D {
    double _Lx;
    double _Ly;
    double _Lz;
public:
    Function3D(double Lx, double Ly, double Lz);
    double a_t();
    double operator()(double x, double y, double z, double t);
};