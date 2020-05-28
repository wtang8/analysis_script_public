#!/bin/python
import numpy as np
import MDAnalysis as mda
import os

folder = 'data/analysis/rg/'

'''
os.system('mkdir %s'%folder)

top = 'prot.tpr'
traj = 'data/xtc_prot/prot_pbc.xtc'
u = mda.Universe(top,traj)
time = np.zeros(len(u.trajectory),dtype=float)
dataA = np.zeros(len(u.trajectory),dtype=float)
dataB = np.zeros(len(u.trajectory),dtype=float)
chainA = u.residues[:46].atoms
chainB = u.residues[46:].atoms

t = 0
for ts in u.trajectory:
    time[t] = u.trajectory.time
    dataA[t] = chainA.atoms.radius_of_gyration()
    dataB[t] = chainB.atoms.radius_of_gyration()
    t += 1

np.savetxt(folder+'rg_time.dat',zip(time,dataA,dataB),fmt=['%.6f','%.6f','%.6f'])
'''

timeh, ht = np.loadtxt('data/hills/HILLS_PTWTE_BIAS.0',usecols=(0,1),unpack=True)
timeh = timeh.astype(int)

timeh = timeh[19::20]
ht = ht[19::20]

cutoff = 25

time, dataA, dataB = np.loadtxt(folder+'rg_time.dat',unpack=True)
time = time.astype(int)
time = time[1:]
dataA = dataA[1:]
dataB = dataB[1:]
wt = np.loadtxt('data/sum_hills/boltzmannweights_int5.dat',usecols=1)

selected = ht < cutoff
time = time[selected]
dataA = dataA[selected]
dataB = dataB[selected]
wt = wt[selected]

print wt.size

def rg_data(data,outfile):
    nblock = 2
    blocksize = time.size/nblock    

    nbin = 41
    bins = np.linspace(0,40,nbin)
    histblock = np.zeros((nblock,nbin-1),dtype=float)

    for iblock in range(nblock):
        wtblock = wt[iblock*blocksize:(iblock+1)*blocksize]
        histblock[iblock,:], bin_edges = np.histogram(data[iblock*blocksize:(iblock+1)*blocksize],bins=bins,normed=True,weights=wtblock)

    hist = np.mean(histblock,axis=0)
    err = np.std(histblock,axis=0)/np.sqrt(nblock)
    np.savetxt(outfile,zip(bins[:-1]+0.5,hist,err),fmt=['%.1f','%.6f','%.6f'])

rg_data(dataA,folder+'rg_hist_A_sub25.dat')
rg_data(dataB,folder+'rg_hist_B_sub25.dat')

