#!/bin/python
import os, math
from numba import jit
import numpy as np, sys
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.lib import mdamath
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

nnode = 16
outfolder = 'corr'
top = 'prot.gro'

dt = 5 # #frames to skip
times = np.arange(0,50001,dt)

#@jit(nopython=True)
def calc_P2(v,nframe):
    nvec = len(v)
    # correlate each component indipendently
    acorr = np.array([np.correlate(v[:,i],v[:,i],'full') for i in range(3)])[:,nvec-1:]
    # sum the correlations for each component
    acorr = np.sum(acorr, axis = 0)
    # divide by the number of values actually measured and return
    acorr /= (nvec - np.arange(nvec))
    # calculate P2 = 0.5*(3*cos(theta)^2-1)
    P2 = np.multiply(acorr,acorr)
    P2 *= 3.
    P2 -= 1.
    P2 /= 2.
    return P2

#def fitfun(x,a,b,c):
#    return a*np.exp(-x/b)+(1-a)*np.exp(-x/c)

def fitfund(x,a,b,c,d):
    return (1-d)*(a*np.exp(-x/b)+(1-a)*np.exp(-x/c))+d

# Parameter for 600 MHz relaxation
#wh = 600.130e6*2*math.pi
#wn = 60.834e6*2*math.pi 

# Parameter for 850 MHz NMR relaxation
wh = 850.13e6*2*math.pi #in Hz
wn = 86.176e6*2*math.pi #in Hz

whmwn = wh-wn #0.764e9*2*math.pi #in Hz
whpwn = wh+wn #0.936e9*2*math.pi #in Hz

hbar = 1.055e-34
planks = 1.055e-34  	#Reduced Planck's constant
permfs = (4.0e-7)*math.pi #Permeability of free space
gammah = 2.6752e8  	#gyromagnetic ratio 1H
gamman = -2.713e7  	#gyromagnetic ratio 15N
gammac = 6.728e7 	#gyromagnetic ratio 13C 

global rnh, rcah, delta_sigma_nh, delta_sigma_cah
rnh = 1.04e-10 	 	#Length of NH bond (Yao et al., JACS, 2010)
rcah = 1.10e-10 		#Length of Ha-13Ca bond (Yamazaki, JACS, 1994)
delta_sigma_nh = -163e-6  #-163 ppm (Yao et al., JACS, 2010)
delta_sigma_cah = 25e-6   #25 ppm (Wei, JACS, 2001)

kappa1 = 0.1*(hbar*permfs*gamman*gammah/2./math.pi/rnh**3)**2
kappa2 = kappa1/2
kappa3 = gammah/gamman*kappa1
k12 = wn*wn*delta_sigma_nh*delta_sigma_nh*(2./15.)
k22 = delta_sigma_nh**2*wn**2/45.

def jomega_err(popt,pcov,omega):
    a = popt[0]
    #d = popt[3]
    t1 = popt[1]
    t2 = popt[2]
    covat1 = t1*t1*pcov[0,0]+a*a*pcov[1,1]+2*a*t1*pcov[0,1]
    covat2 = t2*t2*pcov[0,0]+a*a*pcov[2,2]+2*a*t2*pcov[0,2]
    x1 = 1 + omega*omega*t1*t1
    x2 = 1 + omega*omega*t2*t2
    covfrac1 = covat1/x1/x1+(a*t1*2*omega*t1/x1/x1)**2*pcov[1,1]    
    covfrac2 = covat1/x2/x2+(a*t2*2*omega*t2/x2/x2)**2*pcov[2,2]
    err = (2./5.)*np.sqrt(covfrac1+covfrac2)
    J = (2./5.)*(a*t1/x1+(1-a)*t2/x2)
    return J, err

def r1_err(k11,k12,j0,j0err,jwh,jwherr,jwn,jwnerr,jwhmwn,jwhmwnerr,jwhpwn,jwhpwnerr):
    r1val = k11*(jwhmwn + 3*jwn + 6*jwhpwn) #+ k12*jwn
    r1err = np.sqrt(k11*k11*(jwhmwnerr**2 + (3*jwnerr)**2 + (6*jwhpwnerr)**2)) #+ (k12*jwnerr)**2)
    return r1val, r1err

def r2_err(k21,k22,j0,j0err,jwh,jwherr,jwn,jwnerr,jwhmwn,jwhmwnerr,jwhpwn,jwhpwnerr):
    r2val = k21*(4*j0 + (jwhmwn) + 3*jwn + 6*jwh + 6*jwhpwn) + k22*(4*j0 + 3*jwn)
    r2err = np.sqrt(k21*k21*((4*j0err)**2 + (jwhmwnerr)**2 + (3*jwnerr)**2 + (6*jwherr)**2 + (6*jwhpwnerr)**2) + k22*k22*((4*j0err)**2 + (3*jwnerr)**2))
    return r2val, r2err

def noe_err(k3,jwhpwn,jwhpwnerr,jwhmwn,jwhmwnerr,r1,r1err):
    noeval = 1 + k3*(6*jwhpwn - jwhmwn)/r1
    noeerr = np.sqrt(k3*k3*((6*jwhpwnerr)**2 + jwhmwnerr**2)+r1err*r1err*(noeval-1)**2)/r1
    return noeval, noeerr

w = np.linspace(0,10,num=100000)     

fout = open('relaxation_param.dat','w')
fout.write("#res j0 r1 r2l r2csa NOE\n")


P2ss = [] 
for inode in range(nnode):
     
    traj = 'data/xtc_pbc/prot_nd%d.xtc'%inode # trjconv -pbc mol

    u = mda.Universe(top,traj)    # MDAnalysis reading the files
    rawN = u.select_atoms('name N')
    H = u.select_atoms('name H')
    Hresid = [a.resid for a in H]
    Nresid = []
    N = []
    for a in rawN:
        if a.resid in Hresid:
            Nresid.append(a.resid)
            N.append(a)

    Nindex = np.hstack([N[i].index for i in range(len(N))])
    N = u.atoms[Nindex]
    nres = len(N)

    Hresid = np.array(Hresid)-1
    Nresid = np.array(Nresid)-1
    nresa = int(np.amax(Nresid))+1
    print(nresa,Hresid,Nresid)

    times = np.array([u.trajectory.time for ts in u.trajectory])
    nframe = len(times)
     
    # Displacement vector subtraction with Periodic Boundary Conditions
    vec = np.zeros((nframe,nres,3),dtype=float)
    @jit(parallel=True)
    def cal_vec(t,Npos,Hpos,boxvec):
        for ires in range(nres):
            vec[t,ires] = Npos[ires] - Hpos[ires]
            for i in range(3):
                dot = np.dot(vec[t,ires],boxvec[i])/np.dot(boxvec[i],boxvec[i])
                if dot > .5:
                    vec[t,ires] -= boxvec[i]
                elif dot < -.5:
                    vec[t,ires] += boxvec[i]
            vec[t,ires] = vec[t,ires]/np.sqrt(np.dot(vec[t,ires],vec[t,ires]))
        return 
    
    t = 0

    for ts in u.trajectory:
        boxvec = mdamath.triclinic_vectors(u.dimensions)
        cal_vec(t,N.positions,H.positions,boxvec)
        t += 1
        sys.stdout.write('Time = %.1f\r'%u.trajectory.time)
        sys.stdout.flush()
    np.save(outfolder+'/vec_nd%d.npy'%inode,vec)
    sys.stdout.write('\n')
    
    u = mda.Universe(top)#,traj)    # MDAnalysis reading the files
    rawN = u.select_atoms('name N')
    H = u.select_atoms('name H')
    Hresid = [a.resid for a in H]
    Nresid = []
    N = []
    for a in rawN:
        if a.resid in Hresid:
            Nresid.append(a.resid)
            N.append(a)

    Nindex = np.hstack([N[i].index for i in range(len(N))])
    N = u.atoms[Nindex]
    nres = len(N)

    Hresid = np.array(Hresid)-1
    Nresid = np.array(Nresid)-1
    nresa = int(np.amax(Nresid))+1
    
    nframe = len(times)
    nres = len(Hresid)
    
    P2s = []    
    vec = np.load(outfolder+'/vec_nd%d.npy'%inode)
    vec = vec[::dt,:,:]
    for ires in range(nres):
        resid = Hresid[ires]
        P2 = calc_P2(vec[:,ires,:],nframe)
        P2s.append(P2)
        sys.stdout.write('calculating P2, residue number = %d\r'%resid)
        sys.stdout.flush()
    P2ss.append(P2s)

P2ss = np.array(P2ss)
avgP2 = np.mean(P2ss,axis=0)
np.save('corr/P2.npy',avgP2)

times = times[:2000]
avgP2 = np.load('corr/P2.npy')[:,:2000]

for ires in range(nres):
    res = Hresid[ires]+1
     
    popt, pcov = curve_fit(fitfund,times,avgP2[ires],p0=(0.3,50000.0,50000.0,0.0),bounds=([0,0,0,-1],[1,1e7,1e7,1]))
    yy = fitfund(times,*popt)
    
    plt.plot(times,avgP2[ires])
    plt.plot(times,yy)
    plt.savefig('corr/res%d_corr_fit_d.png'%res)
    plt.close()
     
    popt[1] /= 1e12 #in ns
    popt[2] /= 1e12 #in ns
    pcov[1,:] /= 1e12
    pcov[:,1] /= 1e12
    pcov[2,:] /= 1e12
    pcov[:,2] /= 1e12 

    j0, j0err = jomega_err(popt,pcov,0)
    jwh, jwherr = jomega_err(popt,pcov,wh)
    jwn, jwnerr = jomega_err(popt,pcov,wn)
    jwhmwn, jwhmwnerr = jomega_err(popt,pcov,whmwn)# in s
    jwhpwn, jwhpwnerr = jomega_err(popt,pcov,whpwn) # in s
    
    r1, r1err = r1_err(kappa1,k12,j0,j0err,jwh,jwherr,jwn,jwnerr,jwhmwn,jwhmwnerr,jwhpwn,jwhpwnerr)
    r2, r2err = r2_err(kappa2,k22,j0,j0err,jwh,jwherr,jwn,jwnerr,jwhmwn,jwhmwnerr,jwhpwn,jwhpwnerr)
    NOE, NOEerr = noe_err(kappa3,jwhpwn,jwhpwnerr,jwhmwn,jwhmwnerr,r1,r1err)

    fout.write("%5d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n" % (res, j0, j0err, r1, r1err, r2, r2err, NOE, NOEerr))


