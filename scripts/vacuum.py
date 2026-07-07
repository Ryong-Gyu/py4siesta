#!/usr/bin/env python

import sys, os
import glob
import matplotlib.pyplot as plt
import numpy as np

here = os.getcwd()

# Conversion parameters
epsil0 = np.float64(1/ (4*np.pi))
ang2bohr = np.float64(1.889725989)

def find_nearest(x, e, value):
    idx = (np.abs(x - value)).argmin()
    return e[idx]

def get_PAV():
    
    sysname = glob.glob('*.RHO')[0]
    sysname = sysname.split('.')[0]
    cmd = 'grep ' + '"Total number of electrons: " stdout.txt > elec'
    os.system(cmd)
    f = open('elec', 'r')
    txt = f.readline()
    nelec = float(txt.split()[-1])
    f.close()
    os.system('rm elec')

    f = open('macroave.in', 'w')
    
    f.write('Siesta            # Which code have you used to get the input data?\n')
    f.write('Potential         # Which is the input data used to compute the band offset?\n')
    f.write('%s                # Name of the file where the input data is stored\n' %sysname)
    f.write('1                 # Number of convolutions required to calculate the macro. ave.\n')
    f.write('0                 # First length for the filter function in macroscopic average\n')
    f.write('0                 # Second length for the filter function in macroscopic average\n')
    f.write('%f                # Total charge\n'%nelec)
    f.write('spline            # Type of interpolation')
    
    f.close()
    
    os.system('macroave macroave.in')
    
    X = []
    E = []
    
    with open('%s.PAV' %sysname) as f:
        for i, l in enumerate(f):
            line = l
            word = line.split()
    
            X.append(float(word[0]))
            E.append(float(word[1]))
    f.close()

    X = np.array(X, dtype = np.float64)
    E = np.array(E, dtype = np.float64)

   # os.system('rm macroave.in')

    return X, E

if __name__ == '__main__':

    x,e = get_PAV()
