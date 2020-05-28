#!/bin/python

import numpy as np, scipy, os
from scipy import optimize

# Find reweighted ensemble that better fit the experimental data
# but without deviating much from original ensemble

# This is done by minimizing G=X^2-TS where
# X^2 = chi square = (simulated value - experimental value) / sigma^2
# S = entropy = wi*log(wi/w0) to encourage the ensemble to stay closer to equal weight
# T = theta = strength of the entropic force

nframe = 151    # number of frames used for ensemble refinement
prefixs = ['CA','CB']   # type of chemical shift data used for refinement

# Read experimental data
expd = np.genfromtxt('data/exp-generic.dat',dtype=None,comments='#')
expv = {}   # Experimental values
expe = {}   # Experimental errors, which is not used here, but can be used in the future if needed
data = {}   # Simulation data
namelistexp = []
namelistsim = []
for i in range(len(expd)):
    name = expd[i][0]
    namelistexp.append(name)
    expv[name] = float(expd[i][1])
    expe[name] = float(expd[i][2])
    data[name] = []

# Read simulation data
for prefix in prefixs:
    for i in range(3,43):
        name = '%s%d'%(prefix,i)    # name of the data point = prefix + residue number e.g. CA10 is dCA chemical shift for residue number 10
        filename = 'data/sim-%s-generic.dat' % name # read the corresponding simulation data
        if not os.path.isfile(filename):
            continue
        data[name] = np.loadtxt('data/sim-%s-generic.dat'%name)
        if len(data[name]) != nframe:
            print 'len(data[%s]) != nframe' % name
            exit()
        print name
        namelistsim.append(name)        

# Find valid data from both simulation and experimental, use them to fit chi square
namelist = list(set(namelistsim) & set(namelistexp))

# EJP -> For better performance, don't do dictionary lookups in the objective function
sim_data = np.array([data[name] for name in namelist])
exp_data = np.array([expv[name] for name in namelist])
w_0 = np.ones(nframe) / nframe

# list of theta 
#thetas = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6]
thetas = [0.8]

for theta in thetas:
    # EJP -> turned this skipping off for testing
    if False and os.path.isfile('weights_entropy_%.1f.dat'%theta):
        print 'skipping theta = %.1f' %theta     # Skipping theta if the output file already exist
        continue
    else:
        print 'running theta = %.1f' %theta
    bnd = [(0.0,1.0)]*(nframe)  # bound w between 0 and 1
    sigma = 0.5                 # uncertainty for chi square
    factor = 1.0 / sigma / sigma / len(namelist)    # prefactor in the reduced chi-square
    fcs = open("x2vS_%.1f.dat"%theta,'w')           # files that records all the chisquare and entropy values along the optimization path
    
    def xsq(x):   # return value of G to optimizer
        # EJP -> use numpy sum here
        w = x / np.sum(x)
        # EJP -> don't compute chi-squared in nested loops
        chisq = factor * np.sum((np.sum(w * sim_data, axis=1) - exp_data) ** 2)
        # EJP -> since lim(w->0) w ln w = 0, use nansum to replace nans with 0
        # EJP -> this will create warnings because of the nans but nansum should
        # EJP -> take care of them and prevent them from getting through
        entropy = -np.nansum(w * np.log(w / w_0))
        total = chisq - theta * entropy     # G = X^2 - TS
        fcs.write('%.6e %.6e %.6e\n' % (total,chisq,entropy))
        return total
    
    # Calling the optimizer
    res = optimize.basinhopping(xsq,w_0,
        T=0.05,
        niter_success=50,
        # EJP -> run longer (was testing this)
        niter=500,
        # EJP -> you should probably set stepsize as appropriate, play with T, and take_step
        # EJP -> I have made it smaller here so that the steps taken are actually valid, but
        # EJP -> they should be prevented from going outside the bounds
        stepsize=0.01,
        # EJP -> display moved to basinhopping level and turned off in local minimizer
        disp=True,
        minimizer_kwargs={'method':'L-BFGS-B','bounds':bnd,'options':{'disp': True}})#,seed=0,constraints={'type':'eq','fun':norm})
    print res

    fcs.close()
    
    xr = res.x
    wr = xr/sum(xr)
    np.savetxt('weights_entropy_%.1f.dat'%theta,wr) # Output the optimized weights
    simv = {}
    for name in namelist:
        simv[name] = 0.0
        for i in range(nframe):
            simv[name] += wr[i] * data[name][i]
        print name, expv[name], simv[name]
