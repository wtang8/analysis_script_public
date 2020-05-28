#!/bin/python
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS, NSResults
#from MDAnalysis.analysis import contacts, distances
from MDAnalysis.lib import distances, mdamath

topfile = sys.argv[1] #'prot.tpr'
trajfile = sys.argv[2] #'data/xtc_prot/prot.xtc'
folder = 'data/analysis/cationpi/'

u = mda.Universe(topfile,trajfile)
pires = u.select_atoms('resname TYR PHE HIS')
catres = u.select_atoms('resname ARG LYS NH2')

npigroup = 0
pigroups = []
piresid = []
piresname = []
for res in pires.residues:
    if res.resname in ("TYR","PHE"):
        pigroup = res.CG + res.CD1 + res.CD2 + res.CE1 + res.CE2 + res.CZ
        pigroups.append(pigroup)
        piresid.append(res.resid)
        piresname.append(res.resname)
        npigroup += 1
    elif res.resname == "HIS":
        pigroup = res.CB + res.CG + res.ND1 + res.CE1 + res.NE2 + res.CD2
        pigroups.append(pigroup)
        piresid.append(res.resid)
        piresname.append(res.resname)
        npigroup += 1

ncat = 0
cations = []
catresid = []
catresname = []
for res in catres.residues:
    if res.resname == "ARG":
        cations.append(res.NH1)
        catresid.append(res.resid)
        catresname.append(res.resname)
        cations.append(res.NH1)
        catresid.append(res.resid)
        catresname.append(res.resname)
        ncat += 2
    elif res.resname == "LYS":
        cations.append(res.NZ)
        catresid.append(res.resid)
        catresname.append(res.resname)
        ncat += 1
    elif res.resname == "NH2":
        cations.append(res.N)
        catresid.append(res.resid)
        catresname.append(res.resname)
        ncat += 1


# Displacement vector subtraction with Periodic Boundary Conditions
def vector(posA,posB,boxvec):
    vec = posA - posB
    for i in range(3):
        dot = np.dot(vec,boxvec[i])/np.dot(boxvec[i],boxvec[i])
        if dot > .5:
            vec -= boxvec[i]
        elif dot < -.5:
            vec += boxvec[i]
    return vec

fout = open(folder+'cationpi.dat','w')
for ts in u.trajectory:
    sys.stdout.write('\rTime = %d'%u.trajectory.time) 
    boxvec = mdamath.triclinic_vectors(u.dimensions)
    for ipi in range(npigroup):
        pigroup = pigroups[ipi]
        picenter = pigroup.centroid(pbc=True)
        for icat in range(ncat):
            cation = cations[icat]
            dist = distances.calc_bonds(coords1=cation.position,coords2=picenter,box=u.dimensions)
            if dist < 7.0:
                vec = vector(cation.position,picenter,boxvec)
                v1 = vector(pigroup.positions[0],pigroup.positions[5],boxvec)
                v2 = vector(pigroup.positions[1],pigroup.positions[4],boxvec)
                v3 = vector(pigroup.positions[2],pigroup.positions[3],boxvec)
                n1 = mdamath.normal(v1,v2)
                n2 = mdamath.normal(v3,v1)
                n3 = mdamath.normal(v3,v2)
                n = (n1 + n2 + n3)
                angle = np.degrees(mdamath.angle(n,vec))
                if angle > 120 or angle < 60:
                    fout.write('%.1f %d %s %d %s %.6f %.6f\n'%(u.trajectory.time,catresid[icat],catresname[icat],piresid[ipi],piresname[ipi],angle,dist))

sys.stdout.write('\n')
fout.close()

