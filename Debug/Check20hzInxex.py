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

FIR=np.loadtxt('/home/dparker/Desktop/FIRtest.txt')
#img='/share/dbp2123/dparker/Code/TestSTC/Sequential/Z_slc_STC.nii.gz'
img='/home/dparker/Desktop/TestPaddedMatrix.nii.gz'
imgHR='/home/dparker/Desktop/TestHiresMatrix.nii.gz'
img=nb.load(img).get_data()
imgHR=nb.load(imgHR).get_data()
start='/share/dbp2123/dparker/Code/TestSTC/Sequential/Z_slc.nii'
BOLD='/share/dbp2123/dparker/Code/TestSTC/SuperiorFrontalL20Hz.txt'
BOLD=np.loadtxt(BOLD)

def norm(data):
    data=data-np.mean(data)
    data=data/(np.amax(data)-np.amin(data))
    return(data)








x,y,z=4,10,5
start=nb.load(start).get_data()
start=start[x,y,z,:]
line=img[x,y,z,:]
#line=line-line[0]
padlen=8
TR=2.0
Hf=20.0
#FIR=norm(FIR)
FIR=FIR-FIR[0]
lenF=len(FIR)
tfir=np.arange(lenF)/20.0
lenT=img.shape[-1]
timg=np.arange(lenT)*2-padlen*2

print "LenF:\t{}".format(lenF)
print lenT
line20hz=np.zeros(int((lenT+1-padlen*2)*(Hf*TR)))
print len(line20hz)

toffset=(padlen)*TR

tspan=lenF/Hf
firStart=tspan/2.0

ntotal=int(np.floor(tspan/TR-1))
print ntotal

skip=int(round(Hf*TR))
cumsum=[]
cumtime=[]


skipfact=40






for t20n in range(0,len(line20hz)):
    t20=t20n/Hf#+toffset
    
    shift=np.mod(t20n-skipfact/2.0,skipfact)
    nlow=int(np.ceil((t20+toffset-firStart-shift/Hf)/TR))+1
    nhigh=int(np.floor((t20+toffset+firStart-shift/Hf)/TR))
    
    ntotal=nhigh-nlow
    
    # if ntotal==7:

    #     
    SampPoints=[]
    SampVals=[]
    linevals=[]
    linepoints=[]
    currentsum=0
    nt=0
    point=int(np.round(lenF-shift))-1
    # Shift Backward
    while point>=0:

        SampPoints.append(point/20.0)
        SampVals.append(FIR[point])
        currentsum+=FIR[point]*line[nhigh-nt]
        linevals.append(line[nhigh-nt])
        linepoints.append(nhigh-nt)
        point=int(np.round(point-skip))
        nt+=1
    
    cumsum.append(currentsum)#/np.sum(SampVals))
    cumtime.append(t20)


print 'len cumsum:\t{}'.format(len(cumsum))
print  'imghr shape:'
print imgHR.shape

cppline=cumsum
xx,yy,zz,tt=imgHR.shape
mxcor=0.0
mxx=-1
mxy=-1
mxz=-1
for ix in range(xx):
    for iy in range(yy):
        for iz in range(zz):
            tempcor=np.corrcoef(cppline,imgHR[ix,iy,iz,:])[0][-1]
            if tempcor>mxcor:
                mxcor=tempcor
                mxx=ix
                mxy=iy
                mxz=iz
            

        
        
    
print 'Max Cor {} at: {} {} {}'.format(mxcor,mxx,mxy,mxz)
sys.exit()
