#!/usr/bin/python
import os, sys, numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from matplotlib import rc, rcParams

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

f3d = np.load('contact_noh.npy')
f3d = f3d[:44,44:]
wt = np.loadtxt('../../sum_hills/boltzmannweights.dat',usecols=1)
wt /= np.sum(wt)
mat = np.dot(f3d,wt)

resA, dataA, errA = np.loadtxt('contact_1d_A.dat',usecols=(0,1,2),unpack=True)
resB, dataB, errB = np.loadtxt('contact_1d_B.dat',usecols=(0,1,2),unpack=True)


xbegin = 11
xend = 54
ybegin = 454
yend = 501

ticks1d = [0,0.5,1.0]

fig = plt.figure()

ax1 = fig.add_axes([0.275,0.3,0.6,0.6])
#cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'blue', 'cyan', 'green', 'yellow','red'])
cmap = mpl.cm.get_cmap('gist_heat_r')
im = ax1.imshow(np.transpose(mat), origin='lower', cmap=cmap, interpolation ='none',vmax=0.12, vmin=0, extent=[xbegin-0.5,xend+0.5,ybegin-0.5,yend+0.5], aspect=float(len(dataA))/(len(dataB)))
ax1.set_xticks(np.arange(int(np.ceil(xbegin/5.0)*5),int(round(xend/5.0)*5),5))
ax1.set_yticks(np.arange(int(np.ceil(ybegin/5.0)*5),int(round(yend/5.0)*5),5))
ax1.tick_params(axis='both',labelsize=0)
ax1.minorticks_on()
ax1.zorder=100

cax = fig.add_axes([0.825, 0.3, 0.04, 0.6])
fig.colorbar(im, cax=cax, orientation='vertical')
cax.set_ylabel('contact propensities',fontsize=12)
cax.tick_params(labelsize=10)

axa = fig.add_axes([0.35,0.1,0.45,0.175])
axa.plot(resA,dataA,'ko-',markersize=3)
axa.errorbar(resA,dataA,errA,0,'k',capsize=2)
axa.set_xticks(np.arange(int(np.ceil(xbegin/5.0)*5),int(round(xend/5.0)*5)+1,5))
axa.set_xlim(xbegin-0.5,xend+0.5)
axa.tick_params(axis='x',which='major',bottom=True,top=True)
axa.tick_params(axis='x',which='minor',bottom=True,top=True)
axa.set_yticks(ticks1d)
axa.tick_params(axis='both',labelsize=10)
axa.set_xlabel('FUS 11-54',fontsize=12)
axa.set_ylabel('<# con. res.>',fontsize=10)
axa.tick_params(axis='y',left=True,right=True,labelright=True,labelleft=False)
axa.tick_params(axis='y',which='minor',left=True,right=True)
axa.minorticks_on()

axb = fig.add_axes([0.15,0.3,0.175,0.6])
axb.plot(dataB,resB,'ko-',markersize=3)
axb.errorbar(dataB,resB,0,errB,'k',capsize=2)
axb.set_xticks(ticks1d)
axb.set_yticks(np.arange(int(np.ceil(ybegin/5.0)*5),int(round(yend/5.0)*5)+1,5))
axb.set_ylim(ybegin-0.5,yend+0.5)
axb.tick_params(axis='y',which='major',left=True,right=True)
axb.tick_params(axis='y',which='minor',left=True,right=True)
axb.tick_params(axis='both',labelsize=10)
axb.set_ylabel('RGG3 454-501',fontsize=12)
axb.set_xlabel('<# con. res.>',fontsize=10)
axb.tick_params(axis='x',top=True,bottom=True,labeltop=True,labelbottom=False)
axb.tick_params(axis='x',which='minor',top=True,bottom=True)
axb.minorticks_on()

plt.savefig('fancy2Dcontactmap.pdf')
plt.show()
plt.close()

