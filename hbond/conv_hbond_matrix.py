#!/bin/python
import numpy as np


backboneAtom = ['C','O','N','H','CA','HA']

nres = 92
begin = 50000.
timestep = 20.
nframe = 7501

bbbb = np.zeros((nres,nres,nframe),dtype=bool)
bbsc = np.zeros((nres,nres,nframe),dtype=bool)
scbb = np.zeros((nres,nres,nframe),dtype=bool)
scsc = np.zeros((nres,nres,nframe),dtype=bool)
whwh = np.zeros((nres,nres,nframe),dtype=bool)

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
        whwh[residD,residA,frame] = True
        whwh[residA,residD,frame] = True
        if atomD in backboneAtom:
            if atomA in backboneAtom:
                bbbb[residD,residA,frame] = True
                bbbb[residA,residD,frame] = True
            else:
                bbsc[residD,residA,frame] = True
                scbb[residA,residD,frame] = True
        else:
            if atomA in backboneAtom:
                scbb[residD,residA,frame] = True
                bbsc[residA,residD,frame] = True
            else:
                scsc[residD,residA,frame] = True
                scsc[residA,residD,frame] = True

np.save('whwh.npy',whwh)
np.save('bbbb.npy',bbbb)
np.save('bbsc.npy',bbsc)
np.save('scbb.npy',scbb)
np.save('scsc.npy',scsc)

