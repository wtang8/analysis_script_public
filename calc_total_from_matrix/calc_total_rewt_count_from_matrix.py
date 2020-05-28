#!/bin/python
import sys, numpy as np

mat = np.load(sys.argv[1])
mat = mat[:44,44:,:]
print mat.shape
nres = 44 #mat.shape[0]/2

wt = np.loadtxt('../../sum_hills/boltzmannweights.dat',usecols=1)

nblock = 2
nframe = mat.shape[2]
blocksize = nframe/nblock
total = np.zeros(nblock,dtype=float)
wtsumbl = np.zeros(nblock,dtype=float)
for ibl in range(nblock):
    matblock = mat[:,:,ibl*blocksize:(ibl+1)*blocksize]
    frametotal = np.sum(matblock,axis=(0,1))
    wtblock = wt[ibl*blocksize:(ibl+1)*blocksize]
    wtsumbl[ibl] = np.sum(wtblock)
    total[ibl] = np.dot(frametotal,wtblock)/wtsumbl[ibl]

sumwt = np.sum(wtsumbl)
mean = np.dot(total,wtsumbl)/sumwt
err = np.sqrt(np.dot((total-mean)**2,wtsumbl)/nblock/sumwt/nblock)
print '%.10f %.10f'%(mean,err)

