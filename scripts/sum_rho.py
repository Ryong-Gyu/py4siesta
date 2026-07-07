#!/usr/bin/env python

import siestaio as io
import grid
import os,glob,sys
import numpy as np
import matplotlib.pyplot as plt

ang2bohr = np.float64(1.889725989)

def sum_rho(cell, grid):

    size = np.shape(grid)

    x = np.zeros(size[3])
    y = np.zeros(size[3])
    count = size[1] * size[2] * size[3]

    print(cell,size)

    vol = abs(np.dot(np.cross(cell[0],cell[1]), cell[2]))
    dvol = vol / count

    sumq = 0

    for isp in range(size[0]):
        for ix in range(size[1]):
            for iy in range(size[2]):
                for iz in range(size[3]):
                    sumq += dvol * grid[isp,ix,iy,iz]

    return sumq

files = sys.argv[1]
here = os.getcwd

cell, mesh, grid = io.readGrid(files)


print(f'max: {np.max(grid)}')
print(f'min: {np.min(grid)}')
print(f'avg: {np.mean(grid)}')

sum_rho = sum_rho(cell, grid)

print('Sum Rho: %17.15e'%sum_rho)
