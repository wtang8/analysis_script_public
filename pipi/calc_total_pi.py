#!/bin/python
import numpy as np

wt = np.loadtxt('../../sum_hills/boltzmannweights_int5.dat',usecols=1)

begin = 50000.
timestep = 20.
nframe = 7501
nresa = 44

cutoff = 0.8
count = np.zeros(nframe,dtype=int)
with open('pipi.dat','r') as f:
    for line in f:
        p = filter(None,line.strip().split(' '))
        resa = int(p[3])
        resb = int(p[4])
        if float(p[7]) > cutoff and ((resa <= nresa and resb > nresa) or (resa > nresa and resb <= nresa)):
            time = int(p[0])
            frame = int((time-begin)/timestep)
            count[frame] += 1

nblock = 2
count = count[1:]
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

