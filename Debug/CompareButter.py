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
import pandas as pd
from scipy import signal


def norm(data):
    data=data-np.mean(data)
    data=data/(np.amax(data)-np.amin(data))
    return(data)


# MyFIR=np.loadtxt('/home/dparker/Desktop/FIRtest.txt')
# 
# fs = 18.0  # Sampling frequency
# # Generate the time vector properly
# t = np.arange(1000) / fs
# signala = np.cos(2.0*np.pi*1.8*t) # with frequency of 100
# pl.plot(t, signala, label='a')
# 
# signalb = np.cos(2.0*np.pi*0.1*t) # frequency 20
# pl.plot(t, signalb, label='b')
# 
# signalc = signala + signalb +np.random.rand(len(signalb))
# pl.plot(t, signalc, label='c')
# 
# fc = .2  # Cut-off frequency of the filter
# w = fc / (fs / 2.0) # Normalize the frequency
# b, a = signal.butter(5, w, 'low')
# output = signal.filtfilt(b, a, signalc)
# pl.plot(t, output, label='filtered')
# 
# 
# 
# #Myoutput=signal.fftconvolve(signalc,MyFIR,'same')
# Myoutput=np.convolve(signalc,MyFIR,'same')
# pl.plot(t,Myoutput,label='MyFIR')
# pl.legend()
# pl.figure()
# pl.plot(MyFIR)
# 
# 
# 
# pl.figure()
# 
# scf=np.abs(np.fft.fftshift(np.fft.fft(signalc)))
# btrf=np.abs(np.fft.fftshift(np.fft.fft(output)))
# mof=np.abs(np.fft.fftshift(np.fft.fft(Myoutput)))
# 
# f=np.fft.fftshift(np.fft.fftfreq(len(signalc))*fs)
# pl.plot(f,scf,label='orig')
# pl.plot(f,btrf,label='butter')
# pl.plot(f,mof,label='my')
# pl.legend()
# 
# 
# 
# pl.show()
























Src='/share/studies/Negative_BOLD/Subjects/P00004973/S0001/CheckerBoard/CheckerBoard_P00004973_S0001.nii'
Cor='/home/dparker/Desktop/FsTest.nii.gz'

x=50
y=40
z=7

SrcNii=nb.load(Src)
SrcData=SrcNii.get_data()

CorNii=nb.load(Cor)
CorData=CorNii.get_data()

SrcLine=SrcData[x,y,z,:]
CorLine=CorData[x,y,z,:]

MyFIR=np.loadtxt('/home/dparker/Desktop/FIRtest.txt')

FIRf=np.fft.fftshift(np.fft.fft(MyFIR))


fc=0.2
fs=18
w=fc / (fs / 2.0)



b,a=signal.butter(5,w,'low')


impulse=np.zeros(1000)
impulse[500]=1
impresp=signal.filtfilt(b,a,impulse)

halfimp=impresp[500:]
mFIR=int(np.floor(len(MyFIR)/2))
halffir=MyFIR[mFIR:]


FF=np.fft.fftshift(np.fft.fftfreq(len(SrcLine)))

pl.plot(norm(halfimp),label='butter')
pl.plot(norm(halffir),label='MyFIR')
pl.figure()
imrespf=np.abs(np.fft.fftshift(np.fft.fft(impresp)))
f=np.fft.fftshift(np.fft.fftfreq(len(imrespf)))*fs
mf=np.fft.fftshift(np.fft.fftfreq(len(FIRf)))*fs

pl.semilogy(f,np.abs(np.fft.fftshift(np.fft.fft(impresp))),label='butter')
pl.semilogy(mf,np.abs(FIRf),label='Mine')
pl.legend()

pl.figure()

TrunkButter=impresp[500-mFIR:500+mFIR+1]
pl.plot(TrunkButter,label='truncated butter')
pl.plot(MyFIR,label='mine')
pl.legend()
pl.figure()

TBf=np.fft.fftshift(np.fft.fft(TrunkButter))
pl.semilogy(mf,np.abs(TBf),label='Truncated Butter')
pl.semilogy(mf,np.abs(FIRf),label='MyFIR')
pl.legend()

pl.figure()

fc=0.2
fs=1
w=fc / (fs / 2.0)



b,a=signal.butter(5,w,'low')

BtrLine = signal.filtfilt(b, a, SrcLine)

pl.plot(SrcLine,label='src')
pl.plot(CorLine,label='filtshift')
pl.plot(BtrLine,label='butter')
pl.legend()
pl.title('FullButter')



BLf=np.fft.fftshift(np.fft.fft(BtrLine))
CLf=np.fft.fftshift(np.fft.fft(CorLine))
SLf=np.fft.fftshift(np.fft.fft(SrcLine))

FIRf=np.fft.fftshift(np.fft.fft(MyFIR))


pl.figure()
pl.plot(FF,np.abs(SLf),label='src')
pl.plot(FF,np.abs(CLf),label='filtershift')
pl.plot(FF,np.abs(BLf),label='butter')
pl.legend()
pl.title('FullButter')

pl.figure()
pl.semilogy(np.abs(SLf),label='src')
pl.semilogy(np.abs(CLf),label='filtershift')
pl.semilogy(np.abs(BLf),label='butter')
pl.legend()
pl.title('FullButter')
pl.figure()

##########################
# Truncated Butter
########################
TrunkButter=TrunkButter[::18]
tkbtrline=np.convolve(SrcLine,TrunkButter,'same')
cutL=22
cutR=340
tkbtrline2=tkbtrline[cutL:cutR]
tbr=np.arange(len(SrcLine))
tbr=tbr[cutL:cutR]
pl.plot(norm(SrcLine),label='src')
pl.plot(norm(CorLine),label='filtshift')
pl.plot(tbr,norm(tkbtrline2),label='butter')
pl.legend()
pl.title('trunkbutter')



BLf=np.fft.fftshift(np.fft.fft(tkbtrline))
CLf=np.fft.fftshift(np.fft.fft(CorLine))
SLf=np.fft.fftshift(np.fft.fft(SrcLine))

FIRf=np.fft.fftshift(np.fft.fft(MyFIR))


pl.figure()
pl.plot(np.abs(SLf),label='src')
pl.plot(np.abs(CLf),label='filtershift')
pl.plot(np.abs(BLf),label='butter')
pl.legend()
pl.title('trunkbutter')

pl.figure()
pl.semilogy(np.abs(SLf),label='src')
pl.semilogy(np.abs(CLf),label='filtershift')
pl.semilogy(np.abs(BLf),label='butter')
pl.legend()
pl.title('trunkbutter')
pl.show()

plot(TrunkButter,label='Truncated Butterworth')

