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

    vol = abs(np.dot(np.cross(cell[0],cell[1]), cell[2]))
    dvol = vol / count

    sumq = 0

    for isp in range(size[0]):
        for ix in range(size[1]):
            for iy in range(size[2]):
                for iz in range(size[3]):
                    sumq += dvol * grid[isp,ix,iy,iz]

    return sumq

files1 = sys.argv[1]
files2 = sys.argv[2]

here = os.getcwd

cell, mesh, grid1 = io.readGrid(files1)
cell, mesh, grid2 = io.readGrid(files2)

s1 = sum_rho(cell, grid1)
s2 = sum_rho(cell, grid2)

io.writeGrid('new.RHO',cell, mesh, s2/s1*grid1)
