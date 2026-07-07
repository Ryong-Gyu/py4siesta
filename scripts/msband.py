#!/usr/bin/env python

import os, sys
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import argparse


def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1

def get_MS_band():

    sysbands = glob.glob('*.bands')[0]
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
    
    #
    #   Get high symmetry point
    #
    
        elif i==specialkline:      # band data
            nspecialk = int(words[0])
        elif i>specialkline:
            
            speck.append(float(words[0]))        
            speclabel.append(words[1][1:-1])
    f.close()
    
    #
    #   Find VBM CBM
    #
    
    vbmloc = []
    cbmloc = []
    vbm = MinE
    cbm = MaxE
    
    for i in range(nbands):
        for j in range(nk):
            if E[i][j] <= E_f and E[i][j] > vbm:
                vbm = E[i][j]
                vbmK = k[j]
    
            if E[i][j] >= E_f and E[i][j] < cbm:
                cbm = E[i][j]
                cbmK = k[j]
    
    bandgap = cbm - vbm
    
    #
    #   Find left, right electrode chemical potential
    #
    
    os.system('rm dat')
    CMD = 'grep '+ "'MSDFT: left, right efs'" + ' stdout.txt | tail -n 1'
    CMD = CMD + ' >> dat'
    os.system(CMD)
    f = open('dat', 'r')
    L = f.readline()
    W = L.split()
    f.close()
#    os.system('rm dat')
    print(L)
    print(W)
    left_chemical = float(W[-2])
    right_chemical = float(W[-1])
    
    
    return k, E, nbands, nspin, speck, speclabel, E_f, left_chemical, right_chemical, bandgap
    

def plot_MS_band(K, E, NBANDS, NSPIN, SPEC, LABEL, LCHEM, RCHEM, EG, EMIN, EMAX, KMIN , KMAX):
    
    ######visualization####
    fig = plt.figure(figsize = (3,4))
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


    EAVG = (LCHEM+RCHEM)/2


    c = ['k','r']
    for isp in range(NSPIN):
        for i in range(NBANDS):
            plt.plot(K, E[i]-EAVG, linewidth=2, color=c[isp])
        
    plt.xlim(KMIN,KMAX)
    plt.ylim(EMIN,EMAX)
    
    spec = []
    label = []
    
    idx = 0
    for i in SPEC:
        if KMIN < i and i < KMAX:
            spec.append(SPEC[idx])
            label.append(LABEL[idx])
        
        idx += 1
    #### X-ticks
    #plt.xticks(SPEC,LABEL,fontsize=16)
    labels = [item.get_text() for item in ax.get_xticklabels()]
    empty_string_labels = ['']*len(labels)
    ax.set_xticklabels(empty_string_labels)
    
    #plt.xticks(spec,label,fontsize=12)
 
    for i in range(len(spec)):
        plt.axvline(x=spec[i], color='k', linestyle='--', linewidth=2)
    
    plt.axhline(y=LCHEM-EAVG, color='r', linestyle='--', linewidth=2)
    plt.axhline(y=RCHEM-EAVG, color='b', linestyle='--', linewidth=2)

    plt.savefig('band.png') 
    plt.close()



def main(emin, emax):

    k, e, nbands, nspin, spec, label, ef, lchem, rchem, eg = get_MS_band()

    kmin = min(k)
    kmax = max(k)

    plot_MS_band(k, e, nbands, nspin, spec, label,lchem, rchem, eg, emin, emax, kmin, kmax)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='SIESTA bandstructure visualization')

    parser.add_argument('--emin', type=float, default=-2,
                         help='Minimum energy for bandstructure (eV)')
    parser.add_argument('--emax', type=float, default=+2,
                         help='Minimum energy for bandstructure (eV)')

    args = parser.parse_args()

    main(args.emin,args.emax)



