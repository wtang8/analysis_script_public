#!/bin/python
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn2
import numpy as np, sys

plt.rcParams['font.size'] = 14
seqa = 'TQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSS'
seqb = 'GPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGG'

rta = np.unique(list(seqa))
rtb = np.unique(list(seqb))
nrta = len(rta)
nrtb = len(rtb)
rtaind = [list(rta).index(a) for a in seqa]
rtbind = [list(rtb).index(b) for b in seqb]
rtna = np.zeros(nrta,dtype=int)
rtnb = np.zeros(nrtb,dtype=int)

lena = len(seqa)
lenb = len(seqb)
for i in range(lena):
    rtna[rtaind[i]] += 1
for j in range(lenb):
    rtnb[rtbind[j]] += 1

def write_label(v,which,row,err,colA,colB):
    try:
        if row[colA] > 0.001:
            v.get_label_by_id(which).set_text(r'%2.1f$\pm$%2.1f'%(row[colA]*100,err[colA]*100))
        else:
            v.get_label_by_id(which).set_text('')
    except:
        return

pAraw = np.loadtxt('pair_count.dat').flatten()
epAraw = np.loadtxt('pair_e_count.dat').flatten()

pairs = []
pA = []
epA = []
k = 0
for i in range(nrta):
    for j in range(nrtb):
        if rtna[i] > 3 and rtnb[j] > 3:
            pairs.append('%s%s'%(rta[i],rtb[j]))
            pA.append(pAraw[k])
            epA.append(epAraw[k])
        k += 1

pA = np.array(pA)
epA = np.array(epA)

ind = np.argsort(pA)
pA = pA[ind]
epA = epA[ind]
pairs = np.array(pairs)[ind]

fig = plt.figure(figsize=(8,6))
plt.bar(range(len(ind)),pA,width=0.8,color='black')
plt.errorbar(range(len(ind)),pA,epA,fmt='none',ecolor='black',capsize=6)

plt.xlabel('FUS 11-54 + RGG3 454-501',fontsize=16)
plt.ylabel('<# vdw counts>',fontsize=16)
plt.xticks(range(len(ind)),pairs,rotation='vertical',fontsize=14)
plt.yticks(np.arange(0,16,2),fontsize=14)
plt.xlim(-0.5,len(pA)-0.5)
plt.ylim(0,14)
plt.tight_layout()
plt.savefig('vdw_count_bar.pdf')
plt.show()
plt.close() 

