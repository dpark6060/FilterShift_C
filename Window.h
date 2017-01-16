#ifndef WINDOW_H
#define WINDOW_H
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include "bessel.h"
#include "miscmaths/miscmaths.h"
#include <fstream>
#include <iterator>
#include <algorithm>


class window {
    int PassZero;
    const float PI;
    float cutoff;
    float SamplingRate;
    float Nyq;
    int N;
    float beta;
    std::vector<float> FIR;
    float StopGain;
    double TranWidth;
    public:
        window();
        window(float co, float sr, float sg, double tw, int pz);
        void kaiserord (float,float);
        float kaiser_atten (int,float);
        void kaiser_beta (float);
        void get_window(std::vector<float>, int);
        void print_info();
        std::vector<float> get_fir();
        
        
    
};
#endif /* WINDOW_H */