#!/usr/bin/env python

import siestaio as io
import grid
import os,glob,sys
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import h5py

ang2bohr = np.float64(1.889725989)

def plot(grid):

    fig = plt.figure(figsize = (5,5))
    ax = fig.add_subplot(111)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    ax.set_xticks([])
    ax.set_yticks([])

    im = ax.imshow(grid,cmap='turbo',vmin=np.min(grid),vmax=np.max(grid),interpolation='bilinear')
    cax = fig.add_axes([ax.get_position().x1+0.01, 
                        ax.get_position().y0,0.02,
                        ax.get_position().height
                       ])
    fig.colorbar(im, cax=cax)

ninput = len(sys.argv)
files = sys.argv[1]
axis = sys.argv[2]
if ninput == 4:
    flip = sys.argv[3]
else:
    flip = 'None'

here = os.getcwd

# 250813 MS Jeong - hdf5 compatible
if h5py.is_hdf5(files):
    with h5py.File(files, "r") as f:
        if 'value' not in f:
            raise KeyError(f"'value' dataset not found in {files}")
        g = f['value'][...]
        grid = g[np.newaxis, ...]
else:
    cell, mesh, grid = io.readGrid(files)

grid1 = grid.sum(axis=0)    # spin sum


if flip =='None':
    pass
elif flip=='x':
    grid1 = np.flip(grid1,0)
elif flip=='y':
    grid1 = np.flip(grid1,1)
elif flip=='z':
    grid1 = np.flip(grid1,2)

if axis=='x':
    grid2 = grid1.sum(axis=0)
elif axis=='y':
    grid2 = grid1.sum(axis=1)
elif axis=='z':
    grid2 = grid1.sum(axis=2)

plot(grid2)   

name = files.split('/')[-1]
print(np.max(grid))
#plt.tight_layout()
plt.savefig(f'2dsurf_{name}_{axis}.png')
plt.close()
            

