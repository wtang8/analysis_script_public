#!/bin/python
import numpy as np, sys, os
import MDAnalysis as mda
from MDAnalysis.analysis import distances

top = 'prot.gro' #sys.argv[1]
nnd = 1
os.system('mkdir data/analysis/intra_dist/ -p')
nblock = 2

for node in range(nnd):
    print 'Node %d'%node
    #traj = 'data/xtc_prot/nd%d_prot_eq.xtc'%node #sys.argv[2]
    traj = 'data/xtc_prot/prot_pbc_eq.xtc'
    #traj = sys.argv[1]

    u = mda.Universe(top,traj)
    CA = u.select_atoms('name CA')
    nCA = len(CA)
    nframe = len(u.trajectory)

    intra = np.zeros((nframe,nCA-1),dtype=float)
    n = np.zeros((nframe,nCA-1),dtype=int)
    chaindist = []
    for i in range(nCA-1,0,-1):
        for j in range(i):
            chaindist.append(j)
    chaindist = np.array(chaindist)
    matsize = len(chaindist)

    t = 0
    for frame in u.trajectory:
        mat = distances.self_distance_array(CA.positions,box=u.dimensions,backend='OpenMP')
        for i in range(matsize):
            intra[t,chaindist[i]] += mat[i]
            n[t,chaindist[i]] += 1
        sys.stdout.write('\rFrame %d'%t)
        t += 1

    intra = np.divide(intra,n)
    blocksize = nframe/nblock
    blockmean = np.zeros((nblock,nCA-1),dtype=float)
    for ib in range(nblock):
        blockdata = intra[ib*blocksize:(ib+1)*blocksize,:]
        blockmean[ib,:] = np.mean(blockdata,axis=0)
        print blockmean[ib,:]
    mean = np.mean(blockmean,axis=0)
    err = np.std(blockmean,axis=0)/np.sqrt(nblock)  
    
    nd = range(1,len(mean)+1)
    np.savetxt('data/analysis/intra_dist/intradist_nd%d_blkerr.dat'%node,zip(nd,mean,err),fmt=['%d','%.6f','%.6f'])
    sys.stdout.write('\n')


