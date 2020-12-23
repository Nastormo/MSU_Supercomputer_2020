#include "function.h"

Function3D::Function3D(double Lx, double Ly, double Lz) 
    :_L({Lx, Ly, Lz})
{
    _a_t = M_PI * sqrt(
        1 / pow(_L[0], 2) +
        1 / pow(_L[1], 2) +
        1 / pow(_L[2], 2)
    ); 
}

double Function3D::operator()(double x, double y, double z, double t) {
    return sin((M_PI / _L[0]) * x) * 
           sin((M_PI / _L[1]) * y) * 
           sin((M_PI / _L[2]) * z) * 
           cos(_a_t * t);
}