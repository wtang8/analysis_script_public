#!/bin/python
import numpy as np, sys
import MDAnalysis as mda

folder = 'data/analysis/rg'
'''
# Radius of gyration over time

top = sys.argv[1] #'prot.gro'
traj = sys.argv[2] #'data/xtc_prot/prot_pbc_eq.xtc'

u = mda.Universe(top,traj)
time = np.zeros(len(u.trajectory),dtype=float)
data = np.zeros(len(u.trajectory),dtype=float)
t = 0
for ts in u.trajectory:
    time[t] = u.trajectory.time
    data[t] = u.atoms.radius_of_gyration()
    t += 1

np.savetxt(folder+'/rg_time.dat',zip(time,data),fmt=['%.6f','%.6f'])
'''
time, data = np.loadtxt(folder+'/rg_time.dat',unpack=True)

'''
# Histogram

nblock = 2
blocksize = time.size/nblock
nbin = 81
bins = np.linspace(0,40,nbin)
histblock = np.zeros((nblock,nbin-1),dtype=float)

for iblock in range(nblock):
    histblock[iblock,:], bin_edges = np.histogram(data[iblock*blocksize:(iblock+1)*blocksize],bins=bins,normed=True)

hist = np.mean(histblock,axis=0)
err = np.std(histblock,axis=0)/np.sqrt(nblock)
np.savetxt('rg_hist.dat',zip(bins[:-1]+0.5,hist,err),fmt=['%.1f','%.6f','%.6f'])
'''

# cummulative Histogram

nblock = 2
blocksize = time.size/nblock
nbin = 81
bins = np.linspace(0,40,nbin)
histblock = np.zeros((nblock,nbin-1),dtype=float)

for iblock in range(nblock):
    histbl, bin_edges = np.histogram(data[iblock*blocksize:(iblock+1)*blocksize],bins=bins,normed=True)
    histblock[iblock,:] = np.cumsum(histbl)
    histblock[iblock,:] /= histblock[iblock,-1]

hist = np.mean(histblock,axis=0)
err = np.std(histblock,axis=0)/np.sqrt(nblock)
np.savetxt(folder+'/rg_cumhist.dat',zip(bins[:-1]+0.5,hist,err),fmt=['%.1f','%.6f','%.6f'])

nblock = 2
blocksize = time.size/nblock
meanbl = np.zeros(nblock,dtype=float)
for iblock in range(nblock):
    meanbl[iblock] = np.mean(data[iblock*blocksize:(iblock+1)*blocksize])
mean = np.mean(meanbl)
err = np.std(meanbl)/np.sqrt(nblock)
print '%.3f %.3f'%(mean, err)



