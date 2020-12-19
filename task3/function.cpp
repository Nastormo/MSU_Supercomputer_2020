#include "function.h"

Function3D::Function3D(double Lx, double Ly, double Lz) 
    :_Lx(Lx), _Ly(Ly), _Lz(Lz) 
{}

double Function3D::a_t() {
    return M_PI * sqrt(
        1 / pow(_Lx, 2) +
        1 / pow(_Ly, 2) +
        1 / pow(_Lz, 2)
    ); 
}

double Function3D::operator()(double x, double y, double z, double t) {
    return sin((M_PI / _Lx) * x) * 
           sin((M_PI / _Ly) * y) * 
           sin((M_PI / _Lz) * z) * 
           cos(a_t() * t);
}