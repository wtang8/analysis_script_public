#!/bin/python
import MDAnalysis as mda, sys
import numpy as np
from MDAnalysis.lib.nsgrid import FastNS, NSResults
from MDAnalysis.lib import mdamath
from MDAnalysis.core.topologyobjects import Angle
from datetime import datetime

now = datetime.now()

traj = 'prot_debug.xtc' #sys.argv[1]
top = 'prot.tpr' #sys.argv[2]
outfile = 'hbond_debug.dat' #sys.argv[3]

u = mda.Universe(top,traj)
nframe = len(u.trajectory)

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


donornames = {
        "BB": ['N'],
        "ALA":[],
        "ARG":['NE','NH1','NH2'],
        "ASN":['ND2'],
        "ASP":[],
        "CYS":['SG'],
        "CYH":[],
        "GLN":['NE2'],
        "GLU":[],
        "GLY":[],
        "HIS":['ND1','NE2'],
        "HSD":['ND1'],
        "HSE":['NE2'],
        "HSP":['ND1','NE2'],
        "ILE":[],
        "LEU":[],
        "LYS":['NZ'],
        "MET":[],
        "PHE":[],
        "PRO":[],
        "SER":['OG'],
        "THR":['OG1'],
        "TRP":['NE1'],
        "TYR":['OH'],
        "VAL":[]
        }

acceptornames = {
        "BB": ['O','OC1','OC2'],
        "ALA":[],
        "ARG":[],
        "ASN":['OD1'],
        "ASP":['OD1','OD2'],
        "CYS":[],
        "CYH":[],
        "GLN":['OE1'],
        "GLU":['OE1','OE2'],
        "GLY":[],
        "HIS":['ND1','NE2'],
        "HSD":['NE2'],
        "HSE":['NE1'],
        "HSP":[],
        "ILE":[],
        "LEU":[],
        "LYS":[],
        "MET":['SD'],
        "PHE":[],
        "PRO":[],
        "SER":['OG'],
        "THR":['OG1'],
        "TRP":[],
        "TYR":['OH'],
        "VAL":[]
        }

protein = u.select_atoms('protein')

donors = mda.core.groups.AtomGroup([],u)
protons = mda.core.groups.AtomGroup([],u)
acceptors = mda.core.groups.AtomGroup([],u)
protondonors = mda.core.groups.AtomGroup([],u)

for a in protein.atoms:
    if a.name in donornames["BB"]:
        proton = u.select_atoms('type H and bonded bynum %d'%(a.ix+1))
        if len(proton) != 0:
            donors += a
            protons += proton
            for i in range(len(proton)):
                protondonors += a
    if a.name in donornames[a.resname]:
        proton = u.select_atoms('(type H or type HO) and bonded bynum %d'%(a.ix+1))
        if len(proton) != 0:
            donors += a
            protons += proton
            for i in range(len(proton)):
                protondonors += a
    if a.name in acceptornames["BB"]:
        acceptors += a
    if a.name in acceptornames[a.resname]:
        acceptors += a

outf = open(outfile,'w')
for ts in u.trajectory:   
    time = int(ts.time)
    sys.stdout.flush()
    sys.stdout.write('%d\r'%ts.time)
    #boxvec = mdamath.triclinic_vectors(u.dimensions)
    nsr = FastNS(cutoff=3.0,coords=acceptors.positions,box=u.dimensions,pbc=True,max_gridsize=5000)
    result = nsr.search(protons.positions)
    nb = result.get_pairs()
    dist = result.get_pair_distances()
    ip = 0
    for i,j in nb:
        angle = Angle([protondonors.ix[i],protons.ix[i],acceptors.ix[j]],u).angle(pbc=True)
        if angle > 120.0:
            outf.write('%d %d %d %s %d %s %s %d %s %.8f %.8f\n'%(time,protons.ix[i],acceptors.ix[j],protons[i].resname,protons[i].resid,protons[i].name,acceptors[j].resname,acceptors[j].resid,acceptors[j].name,dist[ip],angle))
        ip += 1
outf.close()
print(datetime.now()-now)

