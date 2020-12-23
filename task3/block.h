#pragma once
#include <vector>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>

#include "function.h"

class Block {
    std::vector<double> _raw;

    std::vector<double> _shift;
    std::vector<int> _size;
    std::vector<int> _min;

public:
    Block(std::vector<int>& size, std::vector<int>& min, 
        std::vector<double>& shift);

    void printBlock() const;
    void printDiff(Function3D &u, double t) const;
    void saveBlock(std::string &str) const;
    
    double getError(Function3D &u, double t) const;

    std::vector<double>& getData() { return _raw; }

    std::vector<int> getSize() { return _size; }
    int getSizeI() { return _size[0]; }
    int getSizeJ() { return _size[1]; }
    int getSizeK() { return _size[2]; }

    int getMinI() { return _min[0]; }
    int getMinJ() { return _min[1]; }
    int getMinK() { return _min[2]; }

    double getShiftX() { return _shift[0]; }
    double getShiftY() { return _shift[1]; }
    double getShiftZ() { return _shift[2]; }

    double getX(int i) const { return (i + _min[0]) * _shift[0]; }
    double getY(int j) const { return (j + _min[1]) * _shift[1]; }
    double getZ(int k) const { return (k + _min[2]) * _shift[2]; }

    double &operator()(int i, int j, int k);
    double operator()(int i, int j, int k) const;

    double& getElem(int i, int j, int k);
    double getValElem(int i, int j, int k) const;

    std::vector<double> getSlice(int axis, int item) const;
    void setSlice(const std::vector<double>& slice, int axis, int item);

    double lap_h(int i, int j, int k) const;
};