#!/usr/bin/env python

import siestaio as io
import grid
import os,glob,sys
import numpy as np
import matplotlib.pyplot as plt
from numba import jit

ang2bohr = np.float64(1.889725989)


def plot(grid):

    fig = plt.figure(figsize = (5,5))
    ax = fig.add_subplot(111)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    ax.set_xticks([])
    ax.set_yticks([])

    im = ax.imshow(grid,cmap='turbo',vmin=0,vmax=0.5,interpolation='bilinear')
    cax = fig.add_axes([ax.get_position().x1+0.01,
                        ax.get_position().y0,0.02,
                        ax.get_position().height
                       ])
    fig.colorbar(im, cax=cax)

ninput = len(sys.argv)

file1 = sys.argv[1]
file2 = sys.argv[2]
axis = sys.argv[3]
if ninput == 5:
    flip = sys.argv[4]
else:
    flip = 'None'

here = os.getcwd
cell, mesh, grid1 = io.readGrid(file1)
cell, mesh, grid2 = io.readGrid(file2)

grid = np.absolute(grid1-grid2)

grid3 = grid.sum(axis=0)


if flip =='None':
    pass
elif flip=='x':
    grid3 = np.flip(grid3,0)
elif flip=='y':
    grid3 = np.flip(grid3,1)
elif flip=='z':
    grid3 = np.flip(grid3,2)

if axis=='x':
    grid3 = grid3.sum(axis=0)
elif axis=='y':
    grid3 = grid3.sum(axis=1)
elif axis=='z':
    grid3 = grid3.sum(axis=2)


name = file1.split('/')[-1]
plot(grid3)   
plt.savefig(f'{name}_diff.png')
plt.close() 
