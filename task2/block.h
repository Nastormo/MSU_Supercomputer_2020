#pragma once
#include <vector>
#include <math.h>

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

    int getSizeI() { return _sizeI; }
    int getSizeJ() { return _sizeJ; }
    int getSizeK() { return _sizeK; }

    int getMinI() { return _minI; }
    int getMinJ() { return _minJ; }
    int getMinK() { return _minK; }

    double getX(int i) { return (i + _minI) * _shiftX; }
    double getY(int j) { return (j + _minJ) * _shiftY; }
    double getZ(int k) { return (k + _minK) * _shiftZ; }

    double &operator()(int i, int j, int k);
    double operator()(int i, int j, int k) const;

    double& getElem(int i, int j, int k);
    double getValElem(int i, int j, int k) const;

    std::vector<double> getDownX();
    std::vector<double> getUpX();
    std::vector<double> getDownY();
    std::vector<double> getUpY();
    std::vector<double> getDownZ();
    std::vector<double> getUpZ();

    void setDownX(std::vector<double>& downX);
    void setUpX(std::vector<double>& upX);
    void setDownY(std::vector<double>& downY);
    void setUpY(std::vector<double>& upY);
    void setDownZ(std::vector<double>& downZ);
    void setUpZ(std::vector<double>& upZ);

    double lap_h(int i, int j, int k) const;
};