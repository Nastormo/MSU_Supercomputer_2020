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

    double _shiftX;
    double _shiftY;
    double _shiftZ;

    int _sizeI;
    int _sizeJ;
    int _sizeK;

    int _minI;
    int _minJ;
    int _minK;

public:
    Block(int sizeI, int sizeJ, int sizeK, 
        int minI, int minJ, int minK, 
        double shiftX, double shiftY, double shiftZ);

    void printBlock() const;
    void printDiff(Function3D &u, double t) const;
    void saveBlock(std::string &str) const;
    
    double getError(Function3D &u, double t) const;

    int getSizeI() { return _sizeI; }
    int getSizeJ() { return _sizeJ; }
    int getSizeK() { return _sizeK; }

    int getMinI() { return _minI; }
    int getMinJ() { return _minJ; }
    int getMinK() { return _minK; }

    double getX(int i) const { return (i + _minI) * _shiftX; }
    double getY(int j) const { return (j + _minJ) * _shiftY; }
    double getZ(int k) const { return (k + _minK) * _shiftZ; }

    double &operator()(int i, int j, int k);
    double operator()(int i, int j, int k) const;

    double& getElem(int i, int j, int k);
    double getValElem(int i, int j, int k) const;

    std::vector<double> getDownI() const;
    std::vector<double> getUpI() const;
    std::vector<double> getDownJ() const;
    std::vector<double> getUpJ() const;
    std::vector<double> getDownK() const;
    std::vector<double> getUpK() const;

    void setDownI(const std::vector<double>& downI);
    void setUpI(const std::vector<double>& upI);
    void setDownJ(const std::vector<double>& downJ);
    void setUpJ(const std::vector<double>& upJ);
    void setDownK(const std::vector<double>& downK);
    void setUpK(const std::vector<double>& upK);

    double lap_h(int i, int j, int k) const;
};