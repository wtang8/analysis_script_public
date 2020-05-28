#!/bin/python
import sys, numpy as np

wt = np.loadtxt('../sum_hills/boltzmannweights.dat',usecols=1)

filenames = ['contact_noh/contact_noh.npy','contact_noh/contact_noh_count.npy','hbond/whwh_count.npy','hbond/bbbb_count.npy','hbond/bbsc_count.npy','hbond/scbb_count.npy','hbond/scsc_count.npy','hbond/whwh.npy','hbond/bbbb.npy','hbond/bbsc.npy','hbond/scbb.npy','hbond/scsc.npy','pipi/whwh_count.npy','pipi/bbbb_count.npy','pipi/bbsc_count.npy','pipi/scbb_count.npy','pipi/scsc_count.npy','pipi/whwh.npy','pipi/bbbb.npy','pipi/bbsc.npy','pipi/scbb.npy','pipi/scsc.npy','cationpi/count.npy','cationpi/matrix_3d.npy','salt/salt.npy','salt/salt_count.npy']
tags = ['# vdw contact residues','# vdw contact atom pairs','Hydrogen bond counts, whole residues','Hydrogen bond counts BB-BB','Hydrogen bond counts, BB-SC','Hydrogen bond counts, SC-BB','Hydrogen bond counts, SC-SC','Hydrogen bond residues, whole residues','Hydrogen bond residues, BB-BB','Hydrogen bond residues, BB-SC','Hydrogen bond residues, SC-BB','Hydrogen bond residues, SC-SC','sp2/pi counts, whole residues','sp2/pi counts, BB-BB','sp2/pi counts, BB-SC','sp2/pi counts, SC-BB','sp2/pi counts, SC-SC','sp2/pi residues, whole residues','sp2/pi residues, BB-BB','sp2/pi residues, BB-SC','sp2/pi residues, SC-BB','sp2/pi residues, SC-SC','cation-pi counts','cation-pi residues','salt bridge','salt bridge counts']
#print len(filenames), len(tags)

for i in range(len(filenames)):
    filename=filenames[i]
    tag = tags[i]
    mat = np.load(filename)
    mat = mat[:44,44:,:]
    
    nblock = 2
    nframe = mat.shape[2]
    blocksize = nframe/nblock
    total = np.zeros(nblock,dtype=float)
    wtsumbl = np.zeros(nblock,dtype=float)
    for ibl in range(nblock):
        matblock = mat[:,:,ibl*blocksize:(ibl+1)*blocksize]
        frametotal = np.sum(matblock,axis=(0,1))
        wtblock = wt[ibl*blocksize:(ibl+1)*blocksize]
        wtsumbl[ibl] = np.sum(wtblock)
        total[ibl] = np.dot(frametotal,wtblock)/wtsumbl[ibl]

    sumwt = np.sum(wtsumbl)
    mean = np.dot(total,wtsumbl)/sumwt
    err = np.sqrt(np.dot((total-mean)**2,wtsumbl)/nblock/sumwt/nblock)
    #print '%s\t%.6f\t%.6f'%(tag,mean,err)
    sys.stdout.write('%.12f %.12f\n'%(mean,err))

