#!/bin/python
import numpy as np


backboneAtom = ['C','O','N','H','CA','HA']

nres = 92
begin = 50000.
timestep = 20.
nframe = 7501

bbbb = np.zeros((nres,nres,nframe),dtype=int)
bbsc = np.zeros((nres,nres,nframe),dtype=int)
scbb = np.zeros((nres,nres,nframe),dtype=int)
scsc = np.zeros((nres,nres,nframe),dtype=int)
whwh = np.zeros((nres,nres,nframe),dtype=int)

with open('hbond.dat','r') as f:
    for line in f:
        p = line.strip().split(' ')
        time = int(p[0])
        residD = int(p[4])
        residA = int(p[7])
        atomD = p[5]
        atomA = p[8]
        atomindD = int(p[1])
        atomindA = int(p[2])
        resnameD = p[3]
        resnameA = p[6]
        '''
        if atomindD < 6:
            residD = 0
            resnameD = 'ACE'
        if atomindD > 590:
            residD += 2
        elif atomindD > 584:
            residD += 1
            resnameD = 'ACE'
        if atomindA < 6:
            residA = 0
            resnameA = 'ACE'
        if atomindA > 590:
            residA += 2
        elif atomindA > 584:
            residA += 1
            resnameA = 'ACE'
        '''
        frame = int((time-begin)/timestep)
        whwh[residD,residA,frame] += 1
        whwh[residA,residD,frame] += 1
        if atomD in backboneAtom:
            if atomA in backboneAtom:
                bbbb[residD,residA,frame] += 1
                bbbb[residA,residD,frame] += 1
            else:
                bbsc[residD,residA,frame] += 1
                scbb[residA,residD,frame] += 1
        else:
            if atomA in backboneAtom:
                scbb[residD,residA,frame] += 1
                bbsc[residA,residD,frame] += 1
            else:
                scsc[residD,residA,frame] += 1
                scsc[residA,residD,frame] += 1

np.save('whwh_count.npy',whwh)
np.save('bbbb_count.npy',bbbb)
np.save('bbsc_count.npy',bbsc)
np.save('scbb_count.npy',scbb)
np.save('scsc_count.npy',scsc)

