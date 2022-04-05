#!/usr/bin/env python

import glob
import scipy.io as sio
import parbl
import stuffr
import matplotlib.pyplot as plt
import os
import numpy as n
import h5py

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
for fi,f in enumerate(fl):
    a=sio.loadmat(fl[len(fl)-fi-1])
    t0,t1=parbl.determine_t0_24(a)
    print("n_samp %d"%(t0-t_prev))
    z=n.array(a["d_raw"][:,0],dtype=n.complex64)
    zd=stuffr.decimate(z,dec=1000)
    plt.plot(zd.real)
    plt.plot(zd.imag)
    plt.show()
    fvecd=n.fft.fftshift(n.fft.fftfreq(len(zd),d=1000/1e6))
    plt.plot(fvecd,10.0*n.log10(n.abs(n.fft.fftshift(n.fft.fft(zd)))**2.0))
    plt.show()
    for i in range(n_ipp):
        
        plt.plot(zd.real)
        plt.plot(zd.imag)
        plt.show()
        z_in=z[(i*ipp):(i*ipp+ipp)]
        plt.plot(fvec,10.0*n.log10(n.fft.fftshift(n.abs(n.fft.fft(z_in))**2.0)))
        plt.show()

