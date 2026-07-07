import time,sys,os,glob
import numpy as np
import subprocess
from scipy import interpolate
import argparse

'''
 DFT-1/2 ion file generator v1.0

 Developer: Kyuhwan Lee
 Description: Calculate self-energy potential and generate modified ion file for DFT-1/2 calc.
 Usage:
          python gion.py [rcut] [neutral atom ae pot. file] [half-ionized atom ae pot. file] [reference ion file]
 Revised:
          2021.09.02: add comments and fix some notations (Ryong-Gyu Lee)
 Revised:
 	  2024.10.03: without [rcut] calculation for isolated atoms (Kaptan Rajput)
 Revised:
      2025.01.16: generalize the code (Ryong-Gyu Lee)
'''

def read_ion(ion):

    # read reference ion files
    cond = subprocess.check_output(f'grep -A 1 "Vna" {ion} | grep -v "Vna" | cut -d" " -f11', shell=True)   # Cutoff
    cond2 = subprocess.check_output(f'grep -A 1 "Vna" {ion} | grep -v "Vna" | cut -d" " -f2', shell=True)   # npts
    numt = subprocess.check_output(f'grep -n Vna {ion} | cut -d : -f 1', shell=True)                        # "Vna" data line
    numtt = subprocess.check_output(f'grep -n Chlocal {ion} | cut -d : -f 1', shell=True)                   # "Chlocal" data line

    start = int(numt)+1    # "Vna" data starting line
    end = int(numtt)-1     # "Vna" data ending line
    dlen = int(end-start)  # Total data points

    vmax = float(cond)
    cnt = int(cond2)  # number of points for ion
    step = vmax/(cnt-1)

    ionf = open(ion,'r')

    ionfs = []

    for i in ionf.readlines():
        ionfs.append(i)

    ionx = np.zeros(end-start)
    iony = np.zeros(end-start)

    for i in range(start,end):
        xy = ionfs[i]
        xy = xy.split()
        ionx[i-start] = xy[0]
        iony[i-start] = xy[1]

    numionline = len(ionfs)


    return ionfs, ionx, iony, start, end, dlen, numionline



def write_ion(out, ionfs, ionx, iony, start, end, dlen, numionline):


    # ion file generation
    e=open(out,'w')

    for i in range(0,start):
        e.write(str(ionfs[i]))

    for i in range(0,dlen):
        ix = format(ionx[i], ".17f")
        iy = format(iony[i], ".15E")
        e.write('    ')
        e.write(ix)
        e.write('       ')
        e.write(iy)
        e.write('\n')
    for i in range(end,numionline):
        e.write(str(ionfs[i]))
    e.close


if __name__=="__main__":


    parser = argparse.ArgumentParser(description='DFT-alpha SIESTA ion generator')

    parser.add_argument('--ion1', type=str, default='./ION',
                         help='Sepecify original SIESTA .ion file')
    parser.add_argument('--ion2', type=str, default='./ION',
                         help='Sepecify alpha SIESTA .ion file')
    parser.add_argument('--ion3', type=str, default='./ION',
                         help='Sepecify alpha SIESTA .ion file')

    parser.add_argument('--out', type=str, default='./NEW',
                         help='Sepecify output SIESTA .ion file')

    args = parser.parse_args()


    ionfs0, ionx0, iony0, start, end, dlen, numionline = read_ion(args.ion1)
    ionfs1, ionx1, iony1, start, end, dlen, numionline = read_ion(args.ion2)
    ionfs2, ionx2, iony2, start, end, dlen, numionline = read_ion(args.ion3)

    ionx = ionx0
    iony = iony0 + (iony1 - iony0) + (iony2 - iony0)

    write_ion(args.out, ionfs0, ionx, iony, start, end, dlen, numionline)
