#!/usr/bin/env python
import nibabel as nb
import numpy as np
import sys
import os
import glob
import matplotlib.pyplot as pl
import matplotlib as mpl
import seaborn as sb
import scipy.stats as sts
import subprocess as sp
import pandas as pd
import seaborn as sb


baseDir='/usr/local/fsl/5.0.7/src/filtershift/Debug/'
Orig=os.path.join(baseDir,'TestPaddedMatrix.nii.gz')
St=os.path.join(baseDir,'TestPaddedMatrix_st.nii.gz')

x,y,z=6,10,5
Tr=1.0

Orig=nb.load(Orig).get_data()
St=nb.load(St).get_data()

oLine=Orig[x,y,z,:]
sLine=St[x,y,z,:]

f,ax=pl.subplots(3,1,figsize=(6,9))

a=ax[0]
a.plot(oLine,label='Original')
a.plot(sLine,label='STC')

a.legend()
a.set_title('LPF Cutoff:0.21Hz, 0.1 transwidth')
a=ax[1]
a.set_xlim([0,0.5])
f=np.fft.fftfreq(len(sLine))
f=np.fft.fftshift(f)

of=np.fft.fftshift(np.fft.fft(oLine))
sf=np.fft.fftshift(np.fft.fft(sLine))

a.semilogy(f,np.abs(of))
a.semilogy(f,np.abs(sf))

a.set_title('Timeseries frequency (log)')

a=ax[2]
fir=np.loadtxt('/home/dparker/Desktop/FIRtest.txt')
fir=fir[::25]
fsamp=len(fir)
f2=np.fft.fftshift(np.fft.fftfreq(fsamp))
#f2=f2*fsamp
fwin=np.fft.fftshift(np.fft.fft(fir))
a.semilogy(f2,np.abs(fwin))
a.set_title('Filter Frequency Response (log)')
a.set_xlim([0,0.5])
pl.subplots_adjust(top=0.95,bottom=0.07,left=0.12,right=0.9,hspace=0.26,wspace=0.20)


pl.figure()
pl.plot(oLine-np.mean(oLine))
pl.plot(sLine)

pl.show()