#!/bin/python
import numpy as np, os

seqa = 'TQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSS'
seqb = 'GPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGG'

lena = len(list(seqa))
lenb = len(list(seqb))

rta = list(np.unique(list(seqa)))
rtb = list(np.unique(list(seqb)))
nrta = len(rta)
nrtb = len(rtb)

rtaind = [rta.index(a) for a in seqa]
rtbind = [rtb.index(b) for b in seqb]


wt = np.loadtxt('../../sum_hills/boltzmannweights.dat',usecols=1)
wt /= np.sum(wt)

nblock = 2

sub = np.load('contact_noh_count.npy')
#sub = np.add(sub[:44,44:],np.transpose(sub[44:,:44],(1,0,2)))
rtpair = np.zeros((nrta,nrtb),dtype=int)

nframe = 7501
blsize = nframe/nblock

pr_A = np.zeros((sub.shape[0],sub.shape[1],nblock),dtype=float)
p_A = np.zeros((nrta,nrtb),dtype=float)
ep_A = np.zeros((nrta,nrtb),dtype=float)

def calc_mean_err(blval,wtblsum):
    sumwt = np.sum(wtblsum)
    mean = np.dot(blval,wtblsum)/sumwt
    err = np.sqrt(np.dot((blval-mean)**2,wtblsum)/nblock/nblock/sumwt)
    return mean, err

wtblsum = np.zeros(nblock,dtype=float)
wtsum = np.zeros((sub.shape[0],sub.shape[1],nblock),dtype=float)
wtbl = []
for ib in range(nblock):
    wtbli = wt[ib*blsize:(ib+1)*blsize]
    wtblsum[ib] = np.sum(wtbli)
    wtbli /= wtblsum[ib]
    wtbl.append(wtbli)
    wtsum[:,:,ib] = np.dot(sub[:,:,ib*blsize:(ib+1)*blsize],wtbli)

for i in range(lena):
    for j in range(lenb):
        mean, err = calc_mean_err(wtsum[i,j],wtblsum)
        p_A[rtaind[i],rtbind[j]] += mean
        ep_A[rtaind[i],rtbind[j]] += err*err
        rtpair[rtaind[i],rtbind[j]] += 1

p_A = np.divide(p_A,rtpair)
p_A[rtpair == 0] = 0
ep_A = np.sqrt(ep_A)
ep_A = np.divide(ep_A,rtpair)
ep_A[rtpair == 0] = 0

np.savetxt('pair_count.dat',p_A)
np.savetxt('pair_e_count.dat',ep_A)

