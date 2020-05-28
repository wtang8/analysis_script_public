#!/usr/bin/python
import sys, numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from matplotlib import rc
#import pylab


f3d = np.load('contact_noh.npy')[:44,44:,:]
wt = np.loadtxt('../../sum_hills/boltzmannweights.dat',usecols=1)
wt /= np.sum(wt)
mat = np.dot(f3d,wt)

seqa = 'TQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSS'
seqb = 'GPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGG'

rta = list(np.unique(list(seqa)))
rtb = list(np.unique(list(seqb)))
nrta = len(rta)
nrtb = len(rtb)

rtinda = [rta.index(a) for a in seqa]
rtindb = [rtb.index(a) for a in seqb]

rtf = np.zeros((nrta,nrtb),dtype=float)
rtn = np.zeros((nrta,nrtb),dtype=float)

rtna = np.zeros(nrta,dtype=float)
rtnb = np.zeros(nrtb,dtype=float)
for i in range(len(seqa)):
    rtna[rtinda[i]] += 1
for j in range(len(seqb)):
    rtnb[rtindb[j]] += 1

'''
for i in range(len(seqa)):
    if rtna[i] < 3:
        mat = np.delete(mat,i,0)
        rta = np.delete(rta,i)
for j in range(len(seqb)):
    if rtnb[j] < 3:
        mat = np.delete(mat,j,1)
        rtb = np.delete(rtb,j)
'''

for i in range(len(seqa)):
    for j in range(len(seqb)):
        rtf[rtinda[i],rtindb[j]] += mat[i,j]
        rtn[rtinda[i],rtindb[j]] += 1

rtpa = np.zeros(nrta,dtype=float)
rtpb = np.zeros(nrtb,dtype=float)

for i in range(nrta):
    for j in range(nrtb):
        rtpa[i] += rtn[i,j]
        rtpb[j] += rtn[i,j]
        rtf[i,j] /= rtn[i,j]

nfr = f3d.shape[2]
nresA = f3d.shape[0]
nresB = f3d.shape[1]
nblk = 2
blksz = nfr/nblk

wtblsum = np.zeros(nblk,dtype=float)
rtfa = np.zeros((nblk,nrta),dtype=float)
rtfb = np.zeros((nblk,nrtb),dtype=float)

for ib in range(nblk):
    wtblsum[ib] = np.sum(wt[ib*blksz:(ib+1)*blksz])
    rwdataA = np.dot(np.sum(f3d[:,:,ib*blksz:(ib+1)*blksz],axis=1),wt[ib*blksz:(ib+1)*blksz])/wtblsum[ib]
    rwdataB = np.dot(np.sum(f3d[:,:,ib*blksz:(ib+1)*blksz],axis=0),wt[ib*blksz:(ib+1)*blksz])/wtblsum[ib]
    for i in range(len(seqa)):
        rtfa[ib,rtinda[i]] += rwdataA[i]      
    for j in range(len(seqb)):
        rtfb[ib,rtindb[j]] += rwdataB[j]
    for i in range(nrta):
        rtfa[ib,i] /= rtpa[i]
    for j in range(nrtb):
        rtfb[ib,j] /= rtpb[j]

sumwt = np.sum(wtblsum)
meanA = np.dot(np.transpose(rtfa),wtblsum)/sumwt
errA = np.zeros(nrta,dtype=float)
for irt in range(nrta):
    errA[irt] = np.sqrt(np.dot((rtfa[:,irt]-meanA[irt])**2,wtblsum)/nblk/nblk/sumwt)

meanB = np.dot(np.transpose(rtfb),wtblsum)/sumwt
errB = np.zeros(nrtb,dtype=float)
for irt in range(nrtb):
    errB[irt] = np.sqrt(np.dot((rtfb[:,irt]-meanB[irt])**2,wtblsum)/nblk/nblk/sumwt)

'''
for i in range(nrta-1,-1,-1):
    if rtna[i] < 3:
        meanA = np.delete(meanA,i)
        errA = np.delete(errA,i)
        rta = np.delete(rta,i)
        rtf = np.delete(rtf,i,0)

for j in range(nrtb-1,-1,-1):
    if rtnb[j] < 3:
        meanB = np.delete(meanB,j)
        errB = np.delete(errB,j)
        rtb = np.delete(rtb,j)
        rtf = np.delete(rtf,j,1)

nrta = len(meanA)
nrtb = len(meanB)
'''
fig = plt.figure()

ax1 = fig.add_axes([0.275,0.3,0.6,0.6])
cmap = mpl.cm.get_cmap('gist_heat_r')
im = ax1.imshow(np.transpose(rtf), origin='lower', cmap=cmap, interpolation ='none',
        vmax=0.016, 
        vmin=0, extent=[-0.5,nrta-0.5,-0.5,nrtb-0.5],aspect=float(nrta)/nrtb)
ax1.set_xticks(range(nrta))
ax1.set_yticks(range(nrtb))
ax1.tick_params(axis='both',labelsize=0)
ax1.zorder=100

cax = fig.add_axes([0.825, 0.3, 0.04, 0.6])
fig.colorbar(im, cax=cax, orientation='vertical')
cax.set_ylabel('<norm. contact propensities>',fontsize=12)
cax.tick_params(labelsize=10)

axa = fig.add_axes([0.35,0.1,0.45,0.175])
axa.bar(range(nrta),meanA,color='k',width=0.4,linewidth=0)
axa.errorbar(x=range(nrta),y=meanA,yerr=errA,ecolor='k',capsize=3,marker=None,fmt='none')
axa.set_xticks(range(nrta))
axa.set_xticklabels(rta)
ticklabelobj = axa.get_xticklabels()
for i in range(nrta):
    if rtna[i] < 3:
        plt.setp(ticklabelobj[i],color='grey')
        axa.bar(i,meanA[i],color='grey',width=0.4,linewidth=0)
        axa.errorbar(x=i,y=meanA[i],yerr=errA[i],ecolor='grey',capsize=3,marker=None,fmt='none')

axa.set_xlim(-0.5,nrta-0.5)
axa.tick_params(axis='x',which='major',bottom=True,top=True)
axa.tick_params(axis='x',which='minor',bottom=False,top=False)

axa.set_yticks([0.0,0.005,0.010])
axa.tick_params(axis='both',labelsize=10)
axa.set_xlabel('FUS 11-54',fontsize=12)
axa.set_ylabel('<con. prop.>',fontsize=10)
axa.tick_params(axis='y',left=True,right=True,labelright=True,labelleft=False)
axa.tick_params(axis='y',which='minor',left=True,right=True)
axa.minorticks_on()

axb = fig.add_axes([0.15,0.3,0.175,0.6])
axb.barh(range(nrtb),meanB,color='k',height=0.4,linewidth=0)
axb.errorbar(x=meanB,y=range(nrtb),xerr=errB,ecolor='k',capsize=3,marker=None,fmt='none')

axb.set_xticks([0.0,0.005,0.01])
axb.set_yticks(range(nrtb))
axb.set_yticklabels(rtb)
ticklabelobj = axb.get_yticklabels()
for j in range(nrtb):
    if rtnb[j] < 3:
        plt.setp(ticklabelobj[j],color='grey')
        axb.barh(j,meanB[j],color='grey',height=0.4,linewidth=0)
        axb.errorbar(x=meanB[j],y=j,xerr=errB[j],ecolor='grey',capsize=3,marker=None,fmt='none')
axb.set_ylim(-0.5,nrtb-0.5)
axb.tick_params(axis='y',which='major',left=True,right=True)
axb.tick_params(axis='y',which='minor',left=False,right=False)

axb.tick_params(axis='both',labelsize=10)
axb.set_ylabel('RGG3 454-501',fontsize=12)
axb.set_xlabel('<con. prop.>',fontsize=10)
axb.tick_params(axis='x',top=True,bottom=True,labeltop=True,labelbottom=False)
axb.tick_params(axis='x',which='minor',top=True,bottom=True)
axb.minorticks_on()

plt.savefig('fancy2drestypenorm.pdf')
plt.show()
plt.close()
