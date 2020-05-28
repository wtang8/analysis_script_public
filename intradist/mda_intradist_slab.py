#!/bin/python
import numpy as np, sys
import MDAnalysis as mda
from MDAnalysis.analysis import distances

top = sys.argv[1]
nnd = 1

for node in range(nnd):
    print 'Node %d'%node
    #traj = 'data/xtc_prot/nd%d_prot_eq.xtc'%node #sys.argv[2]
    traj = sys.argv[2]

    u = mda.Universe(top,traj)
    CA = u.select_atoms('name CA')
    totalnCA = len(CA)
    nchain = 40
    nCA = totalnCA / nchain 
    CA = [CA[nCA*i:nCA*(i+1)] for i in range(nchain)]

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
        for ichain in range(nchain):
            mat = distances.self_distance_array(CA[ichain].positions,box=u.dimensions,backend='OpenMP')
            for i in range(matsize):
                intra[t,chaindist[i]] += mat[i]
                n[t,chaindist[i]] += 1
        sys.stdout.write('\rFrame %d'%t)
        t += 1

    intra = np.divide(intra,n)
    mean = np.mean(intra,axis=0)
    err = np.std(intra,axis=0)#/np.sqrt(intra.shape[0])
    nd = range(1,len(mean)+1)
    np.savetxt(sys.argv[3],zip(nd,mean,err),fmt=['%d','%.6f','%.6f'])



