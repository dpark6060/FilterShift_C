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



def upsample (x, k):
  """
  Upsample the signal to the given ratio using a sinc kernel
  input:
    x   a vector or matrix with a signal in each row
    k   ratio to resample to
    returns
    y   the up or downsampled signal
    when downsampling the signal will be decimated using scipy.signal.decimate
  """

  assert k >= 1, 'k must be equal or greater than 1'

  mn = x.shape
  if len(mn) == 2:
    m = mn[0]
    n = mn[1]
  elif len(mn) == 1:
    m = 1
    n = mn[0]
  else:
    raise ValueError ("x is greater than 2D")

  nn = n * k

  xt = np.linspace (1, n, n)
  xp = np.linspace (1, n, nn)

  return interp (xp, xt, x)

def upsample3 (x, k, workers = None):
  """
  Like upsample, but uses the multi-threaded interp3
  """

  assert k >= 1, 'k must be equal or greater than 1'

  mn = x.shape
  if len(mn) == 2:
    m = mn[0]
    n = mn[1]
  elif len(mn) == 1:
    m = 1
    n = mn[0]
  else:
    raise ValueError ("x is greater than 2D")

  nn = n * k

  xt = np.linspace (1, n, n)
  xp = np.linspace (1, n, nn)

  return interp3 (xp, xt, x, workers)


def interp (xp, xt, x):
  """
  Interpolate the signal to the new points using a sinc kernel
  input:
  xt    time points x is defined on
  x     input signal column vector or matrix, with a signal in each row
  xp    points to evaluate the new signal on
  output:
  y     the interpolated signal at points xp
  """

  mn = x.shape
  if len(mn) == 2:
    m = mn[0]
    n = mn[1]
  elif len(mn) == 1:
    m = 1
    n = mn[0]
  else:
    raise ValueError ("x is greater than 2D")

  nn = len(xp)

  y = np.zeros((m, nn))

  for (pi, p) in enumerate (xp):
    si = np.tile(np.sinc (xt - p), (m, 1))
    y[:, pi] = np.sum(si * x)

  return y.squeeze ()








FIR=np.loadtxt('/home/dparker/Desktop/FIRtest.txt')
#img='/share/dbp2123/dparker/Code/TestSTC/Sequential/Z_slc_STC.nii.gz'
img='/home/dparker/Desktop/TestPaddedMatrix.nii.gz'
img=nb.load(img).get_data()
start='/share/dbp2123/dparker/Code/TestSTC/Sequential/Z_slc.nii'
BOLD='/share/dbp2123/dparker/Code/TestSTC/SuperiorFrontalL20Hz.txt'
BOLD=np.loadtxt(BOLD)
def norm(data):
    data=data-np.mean(data)
    data=data/(np.amax(data)-np.amin(data))
    return(data)



cppline=np.loadtxt('/home/dparker/Desktop/CppTs.txt')

# xx,yy,zz,tt=img.shape
# mxcor=0.0
# mxx=-1
# mxy=-1
# mxz=-1
# for ix in range(xx):
#     for iy in range(yy):
#         for iz in range(zz):
#             tempcor=np.corrcoef(cppline,img[ix,iy,iz,:])[0][-1]
#             if tempcor>mxcor:
#                 mxcor=tempcor
#                 mxx=ix
#                 mxy=iy
#                 mxz=iz
#             
# 
#         
#         
#     
# print 'Max Cor {} at: {} {} {}'.format(mxcor,mxx,mxy,mxz)
# sys.exit()




x,y,z=0,10,5
start=nb.load(start).get_data()
start=start[x,y,z,:]
line=img[x,y,z,:]
#line=line-line[0]
lenFIR=len(FIR)
TR=2.0
Hf=20.0
#FIR=norm(FIR)
FIR=FIR-FIR[0]
padlen=np.ceil(lenFIR/(Hf*TR))

lenF=len(FIR)
tfir=np.arange(lenF)/20.0
lenT=img.shape[-1]
timg=np.arange(lenT)*2-padlen*2

print "LenF:\t{}".format(lenF)
print lenT
line20hz=np.zeros(int((lenT+1-padlen*2)*(Hf*TR)))
print len(line20hz)
#line=norm(line)
#FIR=norm(FIR)
toffset=(padlen)*TR

tspan=lenF/Hf
firStart=tspan/2.0

ntotal=int(np.floor(tspan/TR-1))
print ntotal

skip=int(round(Hf*TR))
cumsum=[]
cumtime=[]
# 
# pl.ion()
# 
# ones=np.ones(3)
# fig=pl.figure()
# ax=fig.add_subplot(111)
# line1=ax.stem(timg,line)
# stemFIR,=ax.plot(ones,ones,'o')
# lineFIR,=ax.plot(ones,ones)
# stemLine,=ax.plot(ones,ones,'*')
# currentline,=ax.plot(ones,ones)
# ax.set_xlim((0,80))

skipfact=40

# for i in range(100):
#     print "np.mod(round(t20n-skipfact/2.0),skipfact):"
#     print "np.mod(round({}),{}):".format(i-skipfact/2.0,skipfact)
#     shift=np.mod(i-skipfact/2.0,skipfact)
#     print shift
# 

# print "TIMESEIRES"
# for i in line:
#     print i
# 
# 
# print "\nFIR"
# 
# 
# for i in FIR:
#     print i






for t20n in range(0,len(line20hz)):
    t20=t20n/Hf
    # 
    # shift=np.mod(t20n+lenF,skipfact)
    # nlow=int(np.ceil((t20+toffset-firStart)/TR))
    # nhigh=int(np.floor((t20+toffset+firStart)/TR))
    #
        
    shift=np.mod(t20n+lenF/2,skipfact)
    nlow=int(np.ceil((t20+toffset-firStart)/TR))
    nhigh=int(np.floor((t20+toffset+firStart)/TR))
    
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
        # print "WHILE LOOP:"
        # print "point:\t{}".format(point)
        # print "FIR[point]:\t{}".format(FIR[point])
        # print "nhigh-nt:\t{}".format(nhigh-nt)
        # print "line[nhigh-nt]:\t{}".format(line[nhigh-nt])
        # print "FIR*line:\t{}".format(FIR[point]*line[nhigh-nt])
        
        SampPoints.append(point/20.0)
        SampVals.append(FIR[point])
        currentsum+=FIR[point]*line[nhigh-nt]
        linevals.append(line[nhigh-nt])
        linepoints.append(nhigh-nt)
        point=int(np.round(point-skip))
        nt+=1
    
    # # out="/home/dparker/Desktop/PyFIR_shift{:.0f}.txt".format(shift)
    # # np.savetxt(out,SampVals)
    # # out="/home/dparker/Desktop/PyFIRt_shift{:.0f}.txt".format(shift)
    # # np.savetxt(out,SampPoints)    
    # # 
    # shift=np.mod(t20n-firStart*20,skipfact)
    # nlow=int(np.ceil((t20+toffset-firStart)/TR))
    # nhigh=int(np.floor((t20+toffset+firStart)/TR))
    # ntotal=nhigh-nlow
    # #Shift Forward
    # point=skip-shift
    # while point<lenF:
    #     
    #     SampPoints.append(point/20.0)
    #     SampVals.append(FIR[point])
    #     currentsum+=FIR[point]*line[nlow+nt]
    #     linevals.append(line[nlow+nt])
    #     linepoints.append(nlow+nt)
    #     point=int(np.round(point+skip))
    #     nt+=1
    # print "Last FIR val: \t{}".format(point+skip)
    # print "shift:\t{}".format(shift)
    # print "t20:\t{}".format(t20)
    # print "t20n:\t{}".format(t20n)
    # print "nlow:\t{}".format(nlow)
    # print "nhigh:\t{}".format(nhigh)
    # print "firStart:\t{}".format(firStart)
    # print "tspan:\t{}".format(tspan)
    # print "ntotal:\t{}".format(ntotal)
    # print "FIRsum:\t{}".format(np.sum(SampVals))
    # print "nt:\t{}".format(nt)
    # # for nt in range(ntotal):
    #     if (nlow+nt<lenT):
    #         firsamp=lenF-(shift+nt*skip)-1
    #         SampPoints.append(firsamp/20.0)
    #         SampVals.append(FIR[firsamp])
    #         currentsum+=FIR[firsamp]*line[nhigh-nt]
    #         linevals.append(line[nlow+nt])
    #         linepoints.append(nlow+nt)
    #         
    
    #cumsum.append(currentsum)
    cumsum.append(currentsum)#/np.sum(SampVals))
    cumtime.append(t20)        
    # if np.mod(t20n,45)==0:
    # #     
    # cumsum.append(currentsum/np.sum(SampVals))
    # cumtime.append(t20)
    # # # #pl.plot(timg,line)
    # # print "Last FIR val: \t{}".format(firsamp)
    # # print "shift:\t{}".format(shift)
    # # print "t20:\t{}".format(t20)
    # # print "t20n:\t{}".format(t20n)
    # # print "nlow:\t{}".format(nlow)
    # # print "nhigh:\t{}".format(nhigh)
    # # print "firStart:\t{}".format(firStart)
    # # print "tspan:\t{}".format(tspan)
    # # print "ntotal:\t{}".format(ntotal)
    # # # print "FIRsum:\t{}".format(np.sum(SampVals))
    # # 
    # stemLine.set_ydata(np.array(linevals))
    # stemLine.set_xdata(np.array(linepoints)*2)
    # 
    # stemFIR.set_ydata(np.array(SampVals))
    # stemFIR.set_xdata(t20+np.array(SampPoints))
    # lineFIR.set_ydata(FIR)
    # lineFIR.set_xdata(t20+tfir)
    # currentline.set_ydata(cumsum)
    # currentline.set_xdata(cumtime)
    
    # pl.stem(t20+np.array(SampPoints),SampVals)
    # pl.plot(t20+tfir,FIR,label='FIR')
    # pl.plot(cumtime,cumsum,label='20hzSummedSamples')
    # pl.xlim((0,80))
    # pl.legend()
    # # pl.show()
    # fig.canvas.draw()
    # fig.canvas.flush_events()
    # #
    #
    
    
    
print "lenF:\t{}".format(lenF)
Py20=np.zeros(len(line)*40)
Py20[::40]=line
Py20=np.convolve(FIR,Py20,'same')[320:]
pl.plot(np.arange(len(Py20))/20.0,norm(Py20),label='PyConvolve')
pl.plot(timg,norm(line),label='Orig')
# # pl.stem(t20+np.array(SampPoints),norm(SampVals),label='lastsamp')
pl.plot(t20+tfir,FIR,label='FIR')
# pl.plot(cumtime,norm(cumsum),label='Final')
pl.plot(np.arange(len(BOLD))*0.05,norm(BOLD),label='Orig')
# pl.plot(np.arange(len(start))*2,start,label='Raw Data')
# 
# #fftResamp=np.fft.ifft(np.fft.fft(start),len(Py20))
# #pl.plot(np.arange(len(Py20))/20.0,norm(np.real(fftResamp)),label='Fft Upsamp')
# sincUpsamp=upsample(start,40)
# pl.plot(np.arange(len(sincUpsamp))*0.05+1.15385,sincUpsamp,label='SincUpsamp')

# pl.xlim((0,80))
pl.plot(np.array(cumtime),norm(cumsum),label='Final')
img=nb.load('/home/dparker/Desktop/TestHiresMatrix.nii.gz').get_data();ts=img[x,y,z,:]
pl.plot(np.array(cumtime),norm(ts),label='cpp')
pl.legend()
print np.sum(FIR)

pl.figure()
N=len(FIR)
pl.semilogy(np.fft.fftshift(np.fft.fftfreq(N)*20),np.abs(np.fft.fftshift(np.fft.fft(FIR))),label='FIR')
N=len(cumsum)
pl.semilogy(np.fft.fftshift(np.fft.fftfreq(N)*20),np.abs(np.fft.fftshift(np.fft.fft(cumsum))),label='Final')
N=len(ts)
pl.semilogy(np.fft.fftshift(np.fft.fftfreq(N)*20),np.abs(np.fft.fftshift(np.fft.fft(ts))),label='Cpp')
N=len(BOLD)
pl.semilogy(np.fft.fftshift(np.fft.fftfreq(N)*20),np.abs(np.fft.fftshift(np.fft.fft(BOLD))),label='Orig')
pl.legend()
pl.show()
