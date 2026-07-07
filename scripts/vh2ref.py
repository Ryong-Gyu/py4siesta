#!/usr/bin/env python

import siestaio as io
import grid
import os,glob,sys
import numpy as np
import matplotlib.pyplot as plt
from numba import jit

ang2bohr = np.float64(1.889725989)
Ry2eV = 13.6056980659


@jit(nopython=True)
def planeaverage(cell, grid):

    size = np.shape(grid)
    x = np.zeros(size[3])
    y = np.zeros(size[3])
    count = size[0] * size[1] * size[2]

    for iz in range(size[3]):

        plane = 0

        for isp in range(size[0]):
            for ix in range(size[1]):
                for iy in range(size[2]):
                    plane += grid[isp,ix,iy,iz]

        x[iz] = iz * (cell[2,2] / mesh[2])
        y[iz] = plane / count

    return x, y

def macroaverage(x, y, l):

    nx = len(x)
    dx = x[1] - x[0]
    nl = int(dx/l)

    ynew = np.zeros(nx)

    for ix in range(len(x)):
        integral = 0
        for il in range(-int(nl/2),int(nl/2),1):
            index = ix + il
            if index < 0:
                index = index + nx
            elif index > nx-1:
                index = index - nx

            integral += y[index]

        ynew[ix] = integral / nl

    return x, ynew

def macro_average(x, y, l):

    dx = x[1] - x[0]  # Grid spacing
    window_size = int(np.round(l / dx))  # Number of points in averaging window

    if window_size % 2 == 0:
        window_size += 1  # Ensure an odd number for symmetric averaging

    half_window = window_size // 2

    # Use periodic boundary conditions for convolution
    y_padded = np.concatenate([y[-half_window:], y, y[:half_window]])

    kernel = np.ones(window_size) / window_size
    y_avg = np.convolve(y_padded, kernel, mode='valid')

    return y_avg


def get_RHO(files):

    rhofile = glob.glob(files)[0]

    cell, mesh, grid = io.readGrid(rhofile)

    nsm = mesh[0] * mesh[1] * mesh[2]
    vol = abs(np.dot(np.cross(cell[0],cell[1]), cell[2]))
    dvol = vol / nsm

    return cell, mesh, grid, dvol

def write_PAV(X, E):

    f = open('PAV.txt','w')

    for i in range(len(X)):
        x = X[i]
        e = E[i]
        f.write(f'{x:17.15f} {e:17.15e}\n')

if __name__=='__main__':

    here = os.getcwd()

    files = sys.argv[1]
    center = float(sys.argv[2])
    cell, mesh, grid = io.readGrid(files)
    x, y = planeaverage(cell, grid)
    center = center * ang2bohr 


    from scipy import interpolate
    f = interpolate.splrep(x, y, s=0)
    ref = interpolate.splev(center, f, der=0)
    print(f'Reference pot: {ref*Ry2eV}')

    ref_vac = interpolate.splev(np.max(x), f, der=0)
    print(f'Vacuum pot: {ref_vac*Ry2eV}')


    write_PAV((x-center)/ang2bohr, y*Ry2eV)
