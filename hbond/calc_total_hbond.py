#!/bin/python
import numpy as np, sys

wt = np.loadtxt('../../sum_hills/boltzmannweights_int5.dat',usecols=1)
print wt.shape

#nresa = 46
begin = 50000.
timestep = 20.
nframe = 7501
natoma = 572

count = np.zeros(nframe,dtype=int)
with open('hbond.dat','r') as f:
    for line in f:
        p = line.strip().split(' ')
        time = int(p[0])
        frame = int((time-begin)/timestep)
        atomindA = int(p[1])
        atomindB = int(p[2])
        if (atomindA > natoma and atomindB < natoma) or (atomindA < natoma and atomindB > natoma):
            count[frame] += 1
with open('hbond_NH2.dat','r') as f:
    for line in f:
        p = line.strip().split(' ')
        time = int(p[0])
        frame = int((time-begin)/timestep)
        atomindA = int(p[1])
        atomindB = int(p[2])
        if (atomindA > natoma and atomindB < natoma) or (atomindA < natoma and atomindB > natoma):
            count[frame] += 1

count = count[1:]
nblock = 2
blocksize = nframe/nblock
total = np.zeros(nblock,dtype=float)
wtsumbl = np.zeros(nblock,dtype=float)
for ibl in range(nblock):
    frametotal = count[ibl*blocksize:(ibl+1)*blocksize]
    wtblock = wt[ibl*blocksize:(ibl+1)*blocksize]
    wtsumbl[ibl] = np.sum(wtblock)
    total[ibl] = np.dot(frametotal,wtblock)/wtsumbl[ibl]

sumwt = np.sum(wtsumbl)
mean = np.dot(total,wtsumbl)/sumwt
err = np.sqrt(np.dot((total-mean)**2,wtsumbl)/nblock/sumwt/nblock)
print '%.2f %.2f'%(mean,err)


