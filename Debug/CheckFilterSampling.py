#!/usr/bin/env python
import nibabel as nb
import numpy as np
import sys
import os
import glob
import matplotlib.pyplot as pl


FIRs=sorted(glob.glob('/home/dparker/Desktop/FIR_shift*'))
shifts=range(len(FIRs))



times=sorted(glob.glob('/home/dparker/Desktop/FIRt_shift*'))

origFIR=np.loadtxt('/home/dparker/Desktop/FIRtest.txt')

tFIR=np.arange(len(origFIR))/20.0
for sh in shifts:
    f='/home/dparker/Desktop/FIR_shift{}.txt'.format(sh)
    pf='/home/dparker/Desktop/PyFIR_shift{}.txt'.format(sh)
    t='/home/dparker/Desktop/FIRt_shift{}.txt'.format(sh)
    pt='/home/dparker/Desktop/PyFIRt_shift{}.txt'.format(sh)
    
    pl.figure()
    fir=np.loadtxt(f)
    time=np.loadtxt(t)
    pfir=np.loadtxt(pf)
    ptime=np.loadtxt(pt)
    
    pl.stem(time,fir,label='Cpp:{}'.format(np.sum(fir)))
    pl.plot(ptime,pfir,'-*',label='py:{}'.format(np.sum(pfir)))
    pl.legend
    pl.plot(tFIR,origFIR)
    
    pl.legend()
    pl.title(sh)
# for f,t in zip(FIRs,times):
#     print f
#     print t
#     pl.figure()
#     fir=np.loadtxt(f)
#     time=np.loadtxt(t)
#     
#     pl.stem(time,fir)
#     pl.plot(tFIR,origFIR)
#     
#     title=os.path.split(f)[-1]
#     title=os.path.splitext(title)[0]
#     pl.title('{}, sum:{}'.format(title,np.sum(fir)))
pl.show()