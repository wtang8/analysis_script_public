#!/bin/python

#    50000   TYR   TYR   3   6  SC  SC    0.71994
with open('pipi.dat','r') as files:
    for line in files:
        p = filter(None,line.strip().split(' '))
        if p[1] == 'TYR' and p[2] == 'TYR' and p[5] == 'SC' and p[6] == 'SC' and float(p[7]) > 0.9:
            print int(p[0])/20-2500, p
    exit()


