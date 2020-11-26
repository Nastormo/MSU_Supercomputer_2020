#include <stdio.h>
#include <math.h>
#include <iostream>


double a_t(double L_x, double L_y, double L_z) {
    return M_PI * sqrt(
        1 / pow(L_x, 2) +
        1 / pow(L_y, 2) +
        1 / pow(L_z, 2)
    ); 
}

double u_analytical(double x, double y, double z, double t, 
                    double L_x, double L_y, double L_z) {
    return sin((M_PI / L_x) * x) * 
           sin((M_PI / L_y) * y) * 
           sin((M_PI / L_z) * z) * 
           cos(a_t(L_x, L_y, L_z) * t);
}

int main(int argc, char** argv) {
    double L_x = std::stod(argv[1]);
    double L_y = std::stod(argv[2]);
    double L_z = std::stod(argv[3]);
    std::cout << "yes" << std::endl;
}