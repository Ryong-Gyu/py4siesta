#!/usr/bin/env python

## LRG   20200217
#################

import os, sys
import matplotlib.pyplot as plt
import numpy as np
import math
import glob

system = glob.glob('*.DM')[0].split('.')[0]
filename = glob.glob('%s.EIG'%system)[0]

#filename = sys.argv[1]

def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1

tlines = file_len(filename)

L = []
W = []

f = open(filename)
for i in range(tlines):
    line = f.readline()
    words = line.split()
    if i == 0: E_f = float(words[0])     # fermi level
    elif i == 1:                         # min and max of E
        neig = int(words[0])
        nspin = int(words[1])
        nkpt = int(words[2])
    L.append(line)
    W.append(words)
f.close()


kblock = int(math.ceil(float(neig*nspin)/10))

E = []
for i in range(nkpt):
    E.append([])
    for j in range(neig*nspin):
        nline = 2+ kblock*i + int(math.floor(float(j)/10))
        row = int(j%10)
        if j<=9:
            E[-1].append(float(W[nline][row+1]))
        else:
            E[-1].append(float(W[nline][row]))

E = np.array(E)


HOMO = []
HOMO1 = []
HOMO2 = []
HOMO3 = []
LUMO = []
LUMO1 = []
LUMO2 = []
LUMO3 = []


for i in range(nkpt):
    for j in range(neig*nspin):
        if E_f <= E[i][j]  and E_f >= E[i][j-1]:
            HOMO.append(E[i][j-1])
            HOMO1.append(E[i][j-2])
            HOMO2.append(E[i][j-3])
            HOMO3.append(E[i][j-4])
            LUMO.append(E[i][j])
            LUMO1.append(E[i][j+1])
            LUMO2.append(E[i][j+2])
            LUMO3.append(E[i][j+3])



homo = max(HOMO)
homo1 = max(HOMO1)
homo2 = max(HOMO2)
homo3 = max(HOMO3)
lumo = min(LUMO)
lumo1 = min(LUMO1)
lumo2 = min(LUMO2)
lumo3 = min(LUMO3)


print("Homo : %f"%homo)
print("Homo1 : %f"%homo1)
print("Homo2 : %f"%homo2)
print("Homo3 : %f"%homo3)

print("LUMO : %f"%lumo)
print("LUMO1 : %f"%lumo1)
print("LUMO2 : %f"%lumo2)
print("LUMO3 : %f"%lumo3)

print("bandgap :")
print(lumo-homo)






