#!/usr/bin/env python

import glob
import scipy.io as sio
import parbl
import stuffr
import matplotlib.pyplot as plt
import os
import numpy as n
import h5py


# sshfs eiscat@goppi.eiscat.uit.no:/ goppi

idir="/home/j/goppi/data/leo_bpark_2.1u_NO@uhf"

# this is how much data is in there

ipp=20000
n_ipp=640
L=n_ipp*ipp
fvec=n.fft.fftshift(n.fft.fftfreq(ipp,d=1/1e6))

fl=glob.glob("%s/20220405*/*.mat"%(idir))
fl.sort()

a=sio.loadmat(fl[0])
# a bit of magic required to figure out start time.
# turns out EISCAT doesn't record exact experiment start time, but it can
# be inferred from approximate data record time and cycle length
t0,t1=parbl.determine_t0_24(a)

ipp_idx=n.arange(ipp)
t_prev=t0
tvec=n.zeros(len(fl))
pvec=n.zeros(len(fl))

T=n.zeros([n_ipp,200],dtype=n.complex64)

for fi,f in enumerate(fl):
    a=sio.loadmat(fl[len(fl)-fi-1])
    t0,t1=parbl.determine_t0_24(a)
    print("n_samp %d"%(t0-t_prev))
    z=n.array(a["d_raw"][:,0],dtype=n.complex64)
    dec=10
    zd=stuffr.decimate(z,dec=dec)
#    plt.plot(zd.real)
 #   plt.plot(zd.imag)
  #  plt.title(stuffr.unix2datestr(t0/1e6))    
   # plt.show()
   # fvecd=n.fft.fftshift(n.fft.fftfreq(len(zd),d=dec/1e6))
   # plt.plot(fvecd,10.0*n.log10(n.abs(n.fft.fftshift(n.fft.fft(zd)))**2.0))
    
   # plt.show()
    ippd=int(ipp/dec)
    tipp=n.zeros(n_ipp)
    for i in range(n_ipp):
        T[i,:]=zd[(i*ippd):(i*ippd+200)]
        tipp[i]=i*ipp*1e-6
        
    plt.plot(tipp,n.unwrap(n.angle(T[:,9])))
    plt.title(stuffr.unix2datestr(t0/1e6))
    plt.grid(markevery=1.0)
    plt.xlabel("Time (s)")
    plt.ylabel("Phase (rad)")
    plt.show()
    plt.pcolormesh(n.angle(T))
    plt.colorbar()
    plt.show()
    

