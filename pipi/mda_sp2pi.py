#!/usr/bin/python
import sys, os, numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import distances
from MDAnalysis.lib.nsgrid import FastNS, NSResults
from MDAnalysis.lib import mdamath

## Determine atom names constituting the sp2 groups in each sidechain and backbone
SC_atoms = {}
SC_atoms['TRP'] = ['CB', 'CG', 'CD1', 'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2']
SC_atoms['TYR'] = ['CB', 'CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']
SC_atoms['PHE'] = ['CB', 'CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']
SC_atoms['HIS'] = ['CB', 'CG', 'ND1', 'CE1', 'NE2', 'CD2']

SC_atoms['ARG'] = ['NE', 'CZ', 'NH1', 'NH2']
SC_atoms['ASP'] = ['CB', 'CG', 'OD1', 'OD2']
SC_atoms['GLU'] = ['CG', 'CD', 'OE1', 'OE2']
SC_atoms['ASN'] = ['CB', 'CG', 'OD1', 'ND2']
SC_atoms['GLN'] = ['CG', 'CD', 'OE1', 'NE2']

BB_atoms = ['CA','C', 'O', 'OC1', 'OC2'] # OC1 and OC2 are for uncapped C-terminal
# N in next residue will be counted below

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

## import trajectory and topology
top = sys.argv[1]
traj = sys.argv[2]
outpfile = sys.argv[3]

u = mda.Universe(top, traj)

## Save sequence and pi-groups into lists
sequence = [u.residues[i].resname for i in range(len(u.residues))]
groups = []
group_res = []
resid = []
loc = []
start = [0]
n_atoms = 0
nres = len(u.residues)
for i, res in enumerate(u.residues):
    if i+1 < nres and u.residues[i+1].resid != 1:
        group = res.atoms.select_atoms('name ' + ' '.join(BB_atoms)) + u.residues[i+1].atoms.select_atoms('name N')
        if len(group) > 0: 
            groups.append(group)
            n_atoms = len(groups[-1].atoms)
            start.append(start[-1]+n_atoms)
            group_res.append(res.resname)
            resid.append(i)
            loc.append('BB')
    if res.resname in SC_atoms.keys():
        group = res.atoms.select_atoms('name ' + ' '.join(SC_atoms[res.resname])) 
        if len(group) > 0:
            groups.append(group)
            n_atoms = len(groups[-1].atoms)
            start.append(start[-1]+n_atoms)
            group_res.append(res.resname)
            resid.append(i)
            loc.append('SC')

start = start[:-1]
end = start[1:]
end.append(start[-1]+n_atoms)

#print groups
pi_index = np.hstack([groups[i].indices for i in range(len(groups))])

## If two groups have at least one atom in contact, calculate the angle between the sp2 planes
outp_str = ''
print(len(groups))

pigroup = u.atoms[pi_index]

ngroup = len(pi_index)
print(len(u.trajectory))

outp = open(outpfile, 'w')

n = np.zeros((len(groups),3),dtype=float)
for ts in u.trajectory:
    boxvec = mdamath.triclinic_vectors(u.dimensions)
    dist = np.zeros((ngroup,ngroup),dtype=bool)
    nsr = FastNS(cutoff=6.0,coords=pigroup.positions,box=u.dimensions,pbc=True,max_gridsize=5000)
    result = nsr.self_search()
    nb = result.get_pairs()
    dist[nb[:,0],nb[:,1]] = True
    for i, grp1 in enumerate(groups):
        v1 = vector(grp1.atoms[0].position,grp1.atoms[1].position,boxvec)
        v2 = vector(grp1.atoms[1].position,grp1.atoms[2].position,boxvec)
        n[i] = np.cross(v1, v2)
        n[i] /= np.linalg.norm(n[i])
    for i, grp1 in enumerate(groups):
        for j, grp2 in enumerate(groups[i+1:]):
            j += i + 1
            contacts = dist[start[i]:end[i],start[j]:end[j]]
            if (resid[j]-resid[i] < 2) and (loc[i] == 'BB') and (loc[j] =='BB'):
                continue  ## Skip backbone-backbone of adjacent residues
            if (contacts.sum() >= 2): # Ensure there are at least two atoms from each group in contact with atom(s) of the other group
                cont_matrix = distances.distance_array(np.vstack([grp1.atoms.positions + n[i]*1.7, grp1.atoms.positions - n[i]*1.7]), np.vstack([grp2.atoms.positions + n[j]*1.7, grp2.atoms.positions - n[j]*1.7]), box=u.dimensions) <= 1.5
                if (cont_matrix.max(0).sum() > 1) and (cont_matrix.max(1).sum() > 1):
                    dot_product = np.abs(np.dot(n[i], n[j]))
                    outp_str = '%8i %5s %5s %3i %3i %3s %3s %10.5f\n' % (ts.time, group_res[i], group_res[j], resid[i], resid[j], loc[i], loc[j], dot_product)
                    outp.write(outp_str)
    
    if ts.time % 1 == 0:
        sys.stdout.write('\r%i'%ts.time)
        sys.stdout.flush()

outp.close()
