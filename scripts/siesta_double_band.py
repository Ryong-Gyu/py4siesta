#!/usr/bin/env python


import os, sys
import matplotlib.pyplot as plt
import numpy as np
import math
import glob


def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1

def get_band(bandfile):

    sysbands = glob.glob(bandfile)[0]
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
    
    #   Get high symmetry point
        elif i==specialkline:      # band data
            nspecialk = int(words[0])
        elif i>specialkline:
            
            speck.append(float(words[0]))        
            speclabel.append(words[1][1:-1])
    f.close()

    #   Find VBM CBM
    vbm_list = []
    cbm_list = []

    for j in range(nk):
        for i in range(nbands*nspin):
            if E[i][j] < E_f and E[i+1][j] > E_f:
                vbm_list.append(E[i][j])
                cbm_list.append(E[i+1][j])
    vbm = max(vbm_list)
    cbm = min(cbm_list)

    
    return k, E, nbands, nspin, speck, speclabel, E_f, vbm, cbm


def plot_band(K, E, NBANDS, NSPIN, EF, C):
    

    ef = EF

    for i in range(NBANDS):
        plt.plot(K, E[i]-ef, linewidth=2, color=C)


if __name__=='__main__':




    # input arguments
    bandfile1 = sys.argv[1] # neutral atom bands file
    bandfile2 = sys.argv[2] # half-ionzied atom bands file

#    emin = float(sys.argv[3])
#    emax = float(sys.argv[4])
    emin = -10
    emax = 10


    # get band structure
    k1, e1, nbands1, nspin1, spec1, label1, ef1, vbm1, cbm1 = get_band(bandfile1)
    k2, e2, nbands2, nspin2, spec2, label2, ef2, vbm2, cbm2 = get_band(bandfile2)

    kmin = min(k1)
    kmax = max(k1)

    # visualization setting
    fig = plt.figure(figsize = (4,3))
    ax = fig.add_subplot(111)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)

    ax.tick_params(axis="x", direction='in', length=6)
    ax.tick_params(axis="y", direction='in', length=6)

    labels = [item.get_text() for item in ax.get_yticklabels()]
    empty_string_labels = ['']*len(labels)
    ax.set_yticklabels(empty_string_labels)

    plt.xlim(kmin,kmax)
    plt.ylim(emin,emax)
    spec = []
    label = []

    idx = 0
    for i in spec1:
        if kmin < i and i < kmax:
            spec.append(spec1[idx])
            label.append(label1[idx])
        idx += 1

    plt.xticks(spec1,label1,fontsize=16)
    labels = [item.get_text() for item in ax.get_xticklabels()]
    empty_string_labels = ['']*len(labels)
    ax.set_xticklabels(empty_string_labels)
    plt.yticks(fontsize=16)

    for i in range(len(spec)):
        plt.axvline(x=spec[i], color='k', linestyle='--', linewidth=1.5)

    # plot band structure
    plot_band(k1, e1, nbands1, nspin1, vbm1, 'k')
    plot_band(k2, e2, nbands2, nspin2, vbm2, 'r')

    plt.tight_layout()
    plt.savefig('band.png')
    plt.close()

