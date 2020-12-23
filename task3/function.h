#pragma once
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <vector> 

class Function3D {
    std::vector<double> _L;
    double _a_t;
public:
    Function3D(double Lx, double Ly, double Lz);
    std::vector<double> getL() const { return _L; }
    double getLx() { return _L[0]; }
    double getLy() { return _L[1]; }
    double getLz() { return _L[2]; }
    double a_t() { return _a_t; }
    double operator()(double x, double y, double z, double t);
};