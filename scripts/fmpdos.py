#!/usr/bin/env python
import glob
import sys, os
import numpy as np
import matplotlib.pylab as plt
import math
import argparse


def pdos(filename, fermi, emin, emax):

    energy = []
    pdos1 = []
    pdos2 = []

    with open(filename) as f:
        for i, l in enumerate(f):
            line = l
            word = line.split()

            if word[0] == '#':
                pass
            else:
                if float(word[0])-fermi > emin and float(word[0])-fermi < emax:
                    if len(word) == 3:
                        energy.append([float(word[0])-fermi])
                        pdos1.append([float(word[1])])
                        pdos2.append([-float(word[2])])
                    else:
                        energy.append([float(word[0])-fermi])
                        pdos1.append([float(word[1])])
                        pdos2.append([0])

    energy = np.array(energy, dtype = float)
    pdos1 = np.array(pdos1, dtype = float)
    pdos2 = np.array(pdos2, dtype = float)

    return energy, pdos1, pdos2



def fermi(emin,emax):

    sysbands = glob.glob('*.EIG')[-1]
    f = open('%s'%sysbands)
    line = f.readline()
    words = line.split()
    E_f = float(words[0]) 

    '''
    DOS = glob.glob('*.DOS')[0]
    Energy, pdos1, pdos2 = pdos(DOS,E_f,emin,emax)

    tol = 1e-10
    dE = 0
    for i in range(len(Energy)):
        if Energy[i] < 0:
            if abs(pdos1[i]) > tol or abs(pdos2[i]) > tol:
                dE = Energy[i]
    '''

    return E_f


def homo(level = 0):

    filename = glob.glob('*.EIG')[-1]

    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    tlines = i+1

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
    LUMO = []

    for i in range(nkpt):
        for j in range(neig*nspin):
            if E_f <= E[i][j]  and E_f >= E[i][j-1]:
                HOMO.append(E[i][j-1-level])
                LUMO.append(E[i][j])

    homo = max(HOMO)
    lumo = min(LUMO)

    return homo




def main(emin, emax, lv):


    tar_list = []
    PDOS = glob.glob('*.PDOS')[0]

    while(1):
        print('\n=====FMPDOS PYTHON ============================================')
        print('    Type orbital indexes for visualization \n')
        print('Usage: (Atomic symbol or index) _ (N quantum #) _ (L quantum #) _ (M quantum #) \n')
        print('If you specify each element by 0, then it will select all')
        word = str(input())
        orbital = word.split('_')
        if word == '0':
            break
        else:
            os.system('rm %s'%word)
            f = open('inp','w')
            f.write('%s\n'%PDOS)
            f.write('%s\n'%word)
            f.write('%s\n'%orbital[0])
            f.write('%s\n'%orbital[1])

            if len(orbital) >= 3:
                f.write('%s\n'%orbital[2])
                if len(orbital) >= 4:
                    f.write('%s'%orbital[3])
            f.close()
            cmd = 'fmpdos < inp'
            os.system(cmd)
            os.system('rm inp')
            tar_list.append(word)


    Ef = fermi(emin, emax)
    homo1 = homo(level = lv)
    homo2 = homo(level = 0)
    print(homo1-homo2)

    DOS = glob.glob('*.DOS')[0]
    Energy, pdos1, pdos2 = pdos(DOS, Ef, emin, emax)

    data = np.hstack((Energy, pdos1))
    data = np.hstack((data, pdos2))

    for itar in range(len(tar_list)):
        e, p1, p2 = pdos(tar_list[itar], Ef, emin, emax)
        data = np.hstack((data,p1))
        data = np.hstack((data,p2))

    np.savetxt('PDOS.csv', data, delimiter=',')



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='SIESTA bandstructure visualization')

    parser.add_argument('--emin', type=float, default=-6,
                         help='Minimum energy for DOS (eV)')
    parser.add_argument('--emax', type=float, default=+2,
                         help='Maximum energy for DOS (eV)')
    parser.add_argument('--level', type=int, default=0,
                         help='The order of reference eigenvalue level')

    args = parser.parse_args()

    main(args.emin,args.emax,args.level)



