#!/bin/python
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS, NSResults

#trajfile = 'prot.xtc'
trajfile = sys.argv[2] #'data/xtc_prot/FUS_120-163_FUS_454-501_amber99sbwsSTQ_ptwte_s1_nd0_eq.xtc'
topfile = sys.argv[1] #'prot.gro'

u = mda.Universe(topfile,trajfile)

#proteinNOH = u.select_atoms('protein and not type H')
backbone = u.select_atoms('(protein or resname NH2) and backbone and not type H')
sidechain = u.select_atoms('(protein or resname NH2) and not backbone and not type H')
bbatomind = backbone.atoms.ix
scatomind = sidechain.atoms.ix
natoms = len(bbatomind) + len(scatomind)
protein = u.select_atoms('(protein or resname NH2) and not type H')
proteinind = protein.atoms.ix
residues = protein.residues
natombb = len(bbatomind)

resid = np.zeros(natoms,dtype=int)
resnum = 0
for a in residues:
    for ind in a.atoms.ix:
        if ind in bbatomind:
            atomind = np.where(proteinind == ind)[0]
            resid[atomind] = resnum
    resnum += 1
for a in residues:
    for ind in a.atoms.ix:
        if ind in scatomind:
            atomind = np.where(proteinind == ind)[0]
            resid[atomind] = resnum
    resnum += 1

#for i in range(len(protein)):
#    print resid[i], protein[i]

frame = 0
nframe = len(u.trajectory)
matrix = np.zeros((resnum,resnum,nframe),dtype=bool)

for ts in u.trajectory:
    sys.stdout.write('\rFrame %d'%frame)
    sys.stdout.flush()
    nsr = FastNS(cutoff=6.0,coords=protein.positions,box=u.dimensions,pbc=True,max_gridsize=20000)
    result = nsr.self_search()
    nb = result.get_pairs()[::2]
    for pair in nb:
        matrix[resid[pair[0]],resid[pair[1]],frame] = True
    frame += 1

np.save(sys.argv[3],matrix)

