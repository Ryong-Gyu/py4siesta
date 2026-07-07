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
    print(size[0])

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
    files2 = sys.argv[2]
    cell, mesh, grid = io.readGrid(files)
    cell2, mesh2, grid2 = io.readGrid(files2)
    x, y = planeaverage(cell, (grid-grid2))
    center = cell[2,2]/2 
    write_PAV((x-center)/ang2bohr, y*ang2bohr**3)
