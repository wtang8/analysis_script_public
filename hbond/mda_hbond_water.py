#!/bin/python
import MDAnalysis as mda, sys
from MDAnalysis.analysis import hbonds

traj = sys.argv[1]
top = sys.argv[2]
outfile = sys.argv[3]

u = mda.Universe(traj,top)

nframe = len(u.trajectory)

start = 0

#for iframe in range(start,nframe):    
#    sys.stdout.flush()
#    sys.stdout.write('%d\r'%iframe)
#    h = []

h = hbonds.HydrogenBondAnalysis(u, 'protein', 'resname SOL', distance = 3.0, angle = 120.0, pbc = True, verbose = True)
h.run(verbose = True)
print h.count_by_time()

