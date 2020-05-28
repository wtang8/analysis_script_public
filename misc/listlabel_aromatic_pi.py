#!/bin/python
import sys

start = [120,372]
seq = ['SSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS',
       'RADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGG']       
#'GPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGG']

## Determine atom names constituting the sp2 groups in each sidechain and backbone
SC_pi = {}
SC_pi['TRP'] = ['CB', 'CG', 'CD1', 'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2']
SC_pi['TYR'] = ['CB', 'CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']
SC_pi['PHE'] = ['CB', 'CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']

sys.stdout.write('[')
ngrouppi = 0
nresa = 634
currentresid = 0
chain = 1
printedlabel = False
with open(sys.argv[1],'r') as f:
    nline = 0
    for line in f:
        nline += 1
        if nline < 3:
            continue
        p = filter(None,line.strip().split(' '))
        atomind = int(p[2])
        if atomind > nresa:
            currentresid = 0
            chain = 2
            nresa = 99999
        resid = int(p[0][:-3])
        if resid > currentresid:
            currentresid = resid
            printedlabel = False
        res = p[0][-3:]
        if res == 'SOL':
            break
        if res in SC_pi.keys():
            if p[1] in SC_pi[res]:
                if not printedlabel:
                    ngrouppi += 1
                    try:
                        sys.stdout.write('\'%d%s\','%(start[chain-1]+resid-1,seq[chain-1][resid-1]))
                    except IndexError:
                        print chain, resid
                    printedlabel = True
sys.stdout.write(']\n')

SC_cation = {}
SC_cation['ARG'] = ['HH11','HH12','HH21','HH22']
SC_cation['LYS'] = ['HZ1','HZ2','HZ3']

sys.stdout.write('[')
ngroupcat = 0
nresa = 634
currentresid = 0
chain = 1
printedlabel = False
with open(sys.argv[1],'r') as f:
    nline = 0
    for line in f:
        nline += 1
        if nline < 3:
            continue
        p = filter(None,line.strip().split(' '))
        atomind = int(p[2])
        if atomind > nresa:
            currentresid = 0
            chain = 2
            nresa = 99999
        resid = int(p[0][:-3])
        if resid > currentresid:
            currentresid = resid
            printedlabel = False
        res = p[0][-3:]
        if res == 'SOL':
            break
        if res in SC_cation.keys():
            if p[1] in SC_cation[res]:
                if not printedlabel:
                    ngroupcat += 1
                    try:
                        sys.stdout.write('\'%d%s\','%(start[chain-1]+resid-1,seq[chain-1][resid-1]))
                    except IndexError:
                        print chain, resid
                    printedlabel = True
sys.stdout.write(']\n')

