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
    const float PI;
    float cutoff;
    int SamplingRate;
    float Nyq;
    int N;
    float beta;
    std::vector<float> FIR;
    float StopGain;
    float TranWidth;
    public:
        window();
        window(float co, int sr, float sg, float tw);
        void kaiserord (float,float);
        float kaiser_atten (int,int);
        void kaiser_beta (float);
        void get_window(std::vector<float>, int);
        void print_info();
        std::vector<float> get_fir();
        
        
    
};
#endif /* WINDOW_H */