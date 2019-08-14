#!/usr/bin/env python
import nibabel as nb
import numpy as np
import sys
import os
import glob
import matplotlib.pyplot as pl
import matplotlib as mpl
import seaborn as sb
import scipy.signal as sp


stopgain=-28
transwidth=0.1
samplingrate=25
PassZero=1
Nyq=samplingrate/2
cutoff=0.21

N,beta=sp.kaiserord(stopgain,transwidth/Nyq)

cutoff1d=[0.0,cutoff/Nyq]

alpha=0.5*(N-1)
m=[]
h=[]
c=[]
FIR=[]

for i in range(N):
    m.append(i-alpha)
    h.append(cutoff1d[1]*np.sinc(cutoff1d[1]*m[i]))
    # If passzero=0, something here

left=cutoff1d[0]
right=cutoff1d[1]
if left==0:
    scalefreq=0.0
elif right==1:
    scalefreq=1
else:
    scalefreq=0.5*(left+right)
s=0
kais=[]
for i in range(N):
    # print sp.bessel(0,beta*np.sqrt(1-np.power((i-alpha)/alpha,2)),analog=True)
    # print beta
    c.append(np.cos(np.pi*m[i]*scalefreq))
    kais.append(np.i0(beta*np.sqrt(1-np.power((i-alpha)/alpha,2)))/np.i0(beta))
    FIR.append(h[i]*kais[i])
    s=s+FIR[i]*c[i]



pl.plot(c,label='c')
pl.plot(h,label='h')
pl.plot(kais,label='kais')
pl.plot(FIR,label='FIR')
pl.legend()
pl.figure()

Fh=np.fft.fftshift(np.fft.fft(h))
Fc=np.fft.fftshift(np.fft.fft(c))
Fkais=np.fft.fftshift(np.fft.fft(kais))
Ffir=np.fft.fftshift(np.fft.fft(FIR))
f=np.fft.fftshift(np.fft.fftfreq(len(Fc))*samplingrate)

pl.semilogy(f,np.abs(Fh),label='Fh')
pl.semilogy(f,np.abs(Fc),label='Fc')
pl.semilogy(f,np.abs(Fkais),label='Fkais')
pl.semilogy(f,np.abs(Ffir),label='Ffir')

pl.legend()
pl.figure()
Cfir=np.loadtxt('/home/dparker/Desktop/FIRtest.txt')
pl.plot(Cfir,label='c++')
pl.plot(FIR,label='Py')
pl.legend()
pl.figure()
pl.plot(Cfir-FIR)
print np.sum(np.abs(Cfir-FIR))
pl.show()
pl.show()














