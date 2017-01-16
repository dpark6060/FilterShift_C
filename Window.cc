#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include "Window.h"
#include "bessel.h"
#include "miscmaths/miscmaths.h"
#include <fstream>
#include <iterator>
#include <algorithm>


//#include "i0.c"

    //CutOffFreq=0.2
    //SamplingRate=Fnew
    //Nyq = SamplingRate/2.0   
    //N, beta = sig.kaiserord(StopGain,TranWidth/Nyq)        
    //BPF = sig.firwin(N, CutOffFreq, window=('kaiser', beta), pass_zero=True, scale=True, nyq=Nyq)
    //

using namespace MISCMATHS;


window::window():PI(3.14159265359){
    PassZero=1;
    cutoff=0.0;
    SamplingRate=0.0;
    StopGain=0.0;
    TranWidth=0.0;
    Nyq=SamplingRate/2.0;
    N=1;
    beta=0;
    FIR.reserve(N);
    
}

window::window(float co, float sr, float sg, double tw, int pz):cutoff(co),SamplingRate(sr),StopGain(sg),TranWidth(tw),PassZero(pz),PI(3.14159265359)
{
    
    Nyq=SamplingRate/2.0;
    N=1;
    beta=0;
    FIR.reserve(N);
    
}

void window::kaiserord (float ripple, float width ) {
//    NiquistRate = SamplingRate/2.0
//	N, beta = sig.kaiserord(StopGain,TranWidth/NiquistRate)
//	taps = sig.firwin(N, CutOffFreq, window=('kaiser', beta), pass_zero=False, scale=True, nyq=NiquistRate)
	//return (taps,NiquistRate)

    
    float A = std::abs(ripple);    
    if ( A < 8 ) {
        std::cout << "Requested maximum ripple attentuation " << A << "is too small for the Kaiser formula." << std::endl;
        beta=0;
    }
    else {
        kaiser_beta(A);
    }
    std::cout<<"Raw N:"<<std::endl;
    std::cout<<(A-7.95)/2.285/(PI*width)+1<<std::endl;
    N= (int) std::ceil((A-7.95)/2.285/(PI*width)+1);
}

void window::kaiser_beta (float a) {
    if ( a > 50 ){
        beta=0.1102 * (a - 8.7);
    }
    else if (a > 21){
        beta = 0.5842 * std::pow((a - 21),0.4) + 0.07886 * (a - 21);
    }
    else {
        beta=0.0;
    }
}

float window::kaiser_atten ( int n, float width ){
    float a;
    a=2.285 * (n - 1) * PI * width + 7.95;
    return a;
}




void window::get_window(std::vector<float> win, int n){
    std::cout << n <<std::endl;
    
    
    
}


void window::print_info(){
    std::cout << "PI=\t"<<PI<<std::endl;
    std::cout << "Cutoff=\t"<<cutoff<<std::endl;
    std::cout << "SR=\t"<< SamplingRate<<std::endl;
    std::cout << "StopGain=\t"<< StopGain<<std::endl;
    std::cout << "TranWidth=\t"<< TranWidth<<std::endl;
    std::cout << "Nyq=\t"<< Nyq<<std::endl;
    std::cout << "N=\t"<<N <<std::endl;
    std::cout << "beta=\t"<<beta <<std::endl;
    std::cout << bessi(0,2) <<std::endl;

    
}

std::vector<float> window::get_fir(){
        //CutOffFreq=0.2
    //SamplingRate=Fnew
    //Nyq = SamplingRate/2.0   
    //N, beta = sig.kaiserord(StopGain,TranWidth/Nyq)        
    //BPF = sig.firwin(N, CutOffFreq, window=('kaiser', beta), pass_zero=True, scale=True, nyq=Nyq)
    //
    std::vector<float> cutoff_1d(2);
    float alpha;
    float left;
    float right;
    float scale_frequency;
    float s;
    //#int pass_nyquist;
    
    
    kaiserord(StopGain,TranWidth/Nyq); // given our StopGain and transition width, calculate the order (Sets N and Beta)
    bool odd=N&1;
    //#pass_nyquist=odd^PassZero;
    
    std::cout<<N<<std::endl;
    std::cout<<odd<<std::endl;
    
    if (PassZero==0&&!odd)
    {
        N++;
        StopGain=kaiser_atten(N,TranWidth/Nyq);
        kaiser_beta(StopGain);
        //kaiserord(StopGain,TranWidth/Nyq);
        print_info();
        std::cout<<"Recalculated N: "<<N<<std::endl;
    }
    
    odd=N&1;
    //#pass_nyquist=odd^PassZero;
    
    std::vector<float> m;
    m.reserve(N);
    std::vector<float> h;
    h.reserve(N);
    std::vector<float> c;
    c.reserve(N);

    
    if (PassZero==1)
    {
        cutoff_1d[0]=0.0;
        cutoff_1d[1]=cutoff/Nyq;
    }
    else
    {
        cutoff_1d[0]=cutoff/Nyq;
        cutoff_1d[1]=1.0;
        std::cout<<"HPF mode"<<std::endl;
    }
    
    
    
    alpha=0.5*(N-1);
    
    for (int i=0;i<N;i++)
    {
        m[i]=(float) i-alpha;
        h[i]=cutoff_1d[1]*sincfn(cutoff_1d[1]*m[i]);
        h[i]=h[(int) i]-cutoff_1d[0]*sincfn(cutoff_1d[0]*m[i]);
    }
    

    
    FIR.resize(N);
    left=cutoff_1d[0];
    right=cutoff_1d[1];
    if (left==0.0)
    {
        scale_frequency = 0.0;
    }
    else if (right == 1)
    {
        scale_frequency = 1.0;
    }
    else
    {
        scale_frequency=0.5*(left+right);
    }
    
    
    
    
    for (int i=0;i<N;i++)
    {
        c[i]=std::cos(M_PI*m[i]*scale_frequency);
        FIR[i]=h[i]*bessi((int) 0,(double) beta*sqrt(1-pow(static_cast<float>((i-alpha)/alpha),2.0)))/bessi(0,beta);
        s+=FIR[i]*c[i];
    }
    
    //for (int i=0;i<N;i++)
    //{
    //    FIR[i]=FIR[i]/s;
    //}
    
    std::ofstream output_file("/home/dparker/Desktop/FIRtest.txt");
    output_file.precision(32);
    std::ostream_iterator<float> output_iterator(output_file,"\n");
    std::copy(FIR.begin(),FIR.end(),output_iterator);

    
    return FIR;
}