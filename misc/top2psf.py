import parmed as pmd, sys

gromacs = pmd.load_file(sys.argv[1],sys.argv[2])
gromacs.save('amber.psf')
#gromacs.save('amber.dcd')


