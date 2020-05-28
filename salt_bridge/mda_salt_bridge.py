#!/bin/python

#SC_cation = {}
#SC_cation['ARG'] = ['HH11','HH12','HH21','HH22']
#SC_cation['LYS'] = ['HZ1','HZ2','HZ3']
#SC_cation['HIS'] = ['H5']

#SC_anion = {}
#SC_anion['ASP'] = ['OD1','OD2'] 
#SC_anion['GLU'] = ['OE1','OE2']


import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS, NSResults
from MDAnalysis.analysis import contacts, distances

topfile = sys.argv[1]
trajfile = sys.argv[2]

u = mda.Universe(topfile,trajfile)
residues = u.residues
cation = u.select_atoms('resname ARG LYS HIS NH2 and name NH1 NH2 NZ N5 N')#HH11 HH12 HH21 HH22 HZ1 HZ2 HZ3 H5')
anion = u.select_atoms('resname ASP GLU and name OD1 OD2 OE1 OE2')

#print cation, anion
#exit()

cresid = np.zeros(len(cation),dtype=int)
aresid = np.zeros(len(anion),dtype=int)

resnum = 0
for a in u.residues:
    for i in a.atoms.ix:
        if i in cation.ix:
            cresid[np.where(cation.ix == i)[0]] = resnum
        if i in anion.ix:
            aresid[np.where(anion.ix == i)[0]] = resnum
    resnum += 1
print cresid, aresid

#nares = np.unique(aresid).size
#ncres = np.unique(cresid).size

frame = 0
nframe = len(u.trajectory)

mat = np.zeros((resnum,resnum,nframe),dtype=bool)
for ts in u.trajectory:
    sys.stdout.write('\rFrame %d'%frame)
    sys.stdout.flush()
    nsr = FastNS(cutoff=6.0,coords=cation.positions,box=u.dimensions,pbc=True,max_gridsize=5000)
    result = nsr.search(anion.positions)
    nb = result.get_pairs()
    for pair in nb:
        mat[aresid[pair[0]],cresid[pair[1]],frame] = True
        mat[cresid[pair[1]],aresid[pair[0]],frame] = True
    frame += 1
sys.stdout.write('\n')

np.save(sys.argv[3],mat)

