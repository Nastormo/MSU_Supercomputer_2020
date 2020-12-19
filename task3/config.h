#pragma once

class Config {
    double _Lx;
    double _Ly;
    double _Lz;
    double _T;
    double _tau;
    int _N;
    int _K;
    int _ndims[3];
    int _argc;
    char** _argv;

public:
    Config(int argc, char** argv);

    void printConfig();

    double getLx()   { return _Lx; };
    double getLy()   { return _Ly; };
    double getLz()   { return _Lz; };
    double getT()    { return _T; };
    double getTau()  { return _tau; };
    int getN()       { return _N; };
    int getK()       { return _K; };
    int getArgc()    { return _argc; };
    int* getNdims()  { return _ndims; }; 
    char** getArgv() { return _argv; }; 

    void setNdims(int* ndims) {
        for (int i = 0; i < 3; i++) {
            _ndims[i] = ndims[i];
        }  
    };   
};