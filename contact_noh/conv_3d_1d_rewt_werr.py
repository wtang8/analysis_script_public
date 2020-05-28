#!/bin/python
import numpy as np, sys

data = np.load(sys.argv[1])

data = data[:44,44:,:]

wt = np.loadtxt('../../sum_hills/boltzmannweights.dat',usecols=1)
wt /= np.sum(wt)

nfr = data.shape[2]
nresA = data.shape[0]
nresB = data.shape[1]
nblk = 2
blksz = nfr/nblk

rwdataA = np.zeros((nblk,nresA),dtype=float)
rwdataB = np.zeros((nblk,nresB),dtype=float)
wtblsum = np.zeros(nblk,dtype=float)
for ib in range(nblk):
    wtblsum[ib] = np.sum(wt[ib*blksz:(ib+1)*blksz])
    rwdataA[ib] = np.dot(np.sum(data[:,:,ib*blksz:(ib+1)*blksz],axis=1),wt[ib*blksz:(ib+1)*blksz])/wtblsum[ib]
    rwdataB[ib] = np.dot(np.sum(data[:,:,ib*blksz:(ib+1)*blksz],axis=0),wt[ib*blksz:(ib+1)*blksz])/wtblsum[ib]

sumwt = np.sum(wtblsum)
meanA = np.dot(np.transpose(rwdataA),wtblsum)/sumwt
errA = np.zeros(nresA,dtype=float)
for ir in range(nresA):
    errA[ir] = np.sqrt(np.dot((rwdataA[:,ir]-meanA[ir])**2,wtblsum)/nblk/nblk/sumwt)

meanB = np.dot(np.transpose(rwdataB),wtblsum)/sumwt
errB = np.zeros(nresB,dtype=float)
for ir in range(nresB):
    errB[ir] = np.sqrt(np.dot((rwdataB[:,ir]-meanB[ir])**2,wtblsum)/nblk/nblk/sumwt)

resA = np.arange(11,55)
resB = np.arange(454,502)
np.savetxt('contact_1d_A.dat',zip(resA,meanA,errA),fmt=["%d","%.6f","%.6f"])
np.savetxt('contact_1d_B.dat',zip(resB,meanB,errB),fmt=["%d","%.6f","%.6f"])

