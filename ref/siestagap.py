#!/usr/bin/env python
import subprocess
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import math
import glob


def get_eigs(filename):

    with open(filename) as f:
        lines = f.readlines()

    ef = float(lines[0].split()[0])
    neig, nspin, nkpt = map(int, lines[1].split())
    kblock = int(math.ceil(float(neig*nspin)/10))

    E = [[float(lines[2 + kblock*i + math.floor(float(j)/10)].split()[1 + j%10 if j < 10 else j%10])
          for j in range(neig * nspin)]
          for i in range(nkpt)
        ]

    E = np.array(E).reshape((nkpt, nspin, neig))

    return E, ef


def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1

def get_vacuum(filename):

    etot = subprocess.run(f"grep 'dhscf: Vacuum level (max, mean) =' {filename} | awk '{{print $(NF-1)}}'",
                          shell=True,
                          capture_output=True,
                          text=True)

    return float(etot.stdout.strip())



def get_bands(filename):

    sysbands = glob.glob(filename)[0]
    tlines = file_len('%s' %sysbands)

    speck = []
    speclabel = []

    f = open('%s'%sysbands)
    for i in range(tlines):
        line = f.readline()
        words = line.split()
        if i == 0: E_f = float(words[0])     # fermi level
        elif 0<i<2: pass
        elif i == 2:                         # min and max of E
            MinE = np.float64(words[0])
            MaxE = np.float64(words[1])
        elif i == 3:
            nbands = int(words[0])           # number of bands
            nspin = int(words[1])            # number of spins
            nk = int(words[2])               # number of k
            specialkline = int(math.ceil(float(nbands*nspin)/10)*nk + 4)  # line where specialK loc is

            k = np.zeros((nk), dtype = np.float64)
            E = np.zeros((nbands*nspin,nk), dtype = np.float64)

        elif 3 < i < specialkline:

            klabel = (i-4)//int(math.ceil(float(nbands*nspin)/10))

            if (i-4)%int(math.ceil(float(nbands*nspin)/10))==0:
                k[(i-4)//int(math.ceil(float(nbands*nspin)/10))] = np.float64(words[0])
                for j in range(1,len(words)):
                    E[j-1][klabel]=np.float64(words[j])

            else:
                blabel = ((i-4)%int(math.ceil(float(nbands*nspin)/10)))*10
                for j in range(len(words)):
                    E[blabel+j][klabel]=np.float64(words[j])

        elif i==specialkline:      # band data
            nspecialk = int(words[0])
        elif i>specialkline:

            speck.append(float(words[0]))
            speclabel.append(words[1][1:-1])
    f.close()


    occ = fermi(E, E_f)
    occ_mask = E <= E_f
    unocc_mask = E > E_f


#    occ_mask = occ >= 0.01
#    unocc_mask = occ < 0.01

    occ_E = E[occ_mask]
    unocc_E = E[unocc_mask]

    vbm = np.max(occ_E)
    cbm = np.min(unocc_E)

    return E_f, vbm, cbm

def fermi(E, ef, kb = 8.617e-5, T = 300):

    return 1 / (1 + np.exp((E-ef)/ (kb * T)) )    


def get_level(E, ef):

    nkpt, nspin, neig = np.shape(E)

    occ = fermi(E, ef)
    occ_mask = occ >= 0.01
    unocc_mask = occ < 0.01

    occ_E = E[occ_mask] 
    unocc_E = E[unocc_mask]

    homo = np.max(occ_E)
    lumo = np.min(unocc_E)

    return homo, lumo

if __name__== '__main__':

    if sys.argv[1].split('.')[-1] == 'EIG':
        e, ef = get_eigs(sys.argv[1])
        homo, lumo = get_level(e,ef)
        gap = lumo - homo
    elif sys.argv[1].split('.')[-1] == 'bands':
        ef, homo, lumo = get_bands(sys.argv[1]) 
        gap = lumo -homo

    print(f'Fermi level: {ef: 17.15f} eV')
    print(f'HOMO level: {homo: 17.15} eV')
    print(f'LUMO level: {lumo: 17.15} eV')
    print(f'Gap: {gap: 17.15} eV')
