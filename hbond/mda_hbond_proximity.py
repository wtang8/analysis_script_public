#!/bin/python
import MDAnalysis as mda, sys
from MDAnalysis.analysis import hbonds

traj = sys.argv[1]
top = sys.argv[2]
outfile = sys.argv[3]

u = mda.Universe(traj,top)
nframe = len(u.trajectory)

start = 0

for iframe in range(start,nframe):    
    sys.stdout.flush()
    sys.stdout.write('%d\r'%iframe)
    h = []
    h = hbonds.HydrogenBondAnalysis(u, 'protein or resname NH2', 'protein or resname NH2', distance = 3.0, angle = 0.0, pbc = True)
#    h = hbonds.HydrogenBondAnalysis(u, 'resname NH2', 'protein or resname NH2', distance = 3.0, angle = 120.0, pbc = True)
    h.run(start = iframe, stop = iframe+1)
    h.generate_table()
    outf = open(outfile,'a')
    for hbond in h.table:
        outf.write('%d %d %d %s %d %s %s %d %s %.8f %.8f\n'%(hbond[0],hbond[1],hbond[2],hbond[3],hbond[4],hbond[5],hbond[6],hbond[7],hbond[8],hbond[9],hbond[10]))
    outf.close()
