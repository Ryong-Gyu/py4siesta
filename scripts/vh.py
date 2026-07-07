#!/usr/bin/env python

import siestaio as io
import grid
import os,glob,sys
import numpy as np
import matplotlib.pyplot as plt

ang2bohr = np.float64(1.889725989)

vacuum = float(sys.argv[2])
fermi = float(sys.argv[3])

Ry2eV = 13.6056980659 


def set_plot():

    fig = plt.figure(figsize = (5,2))
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

    labels = [item.get_text() for item in ax.get_xticklabels()]
    empty_string_labels = ['']*len(labels)
    ax.set_xticklabels(empty_string_labels)

def plot_PAV(X, E, c='r'):


    plt.xlim([X[0],X[-1]])
    plt.ylim([-15,1])
    plt.plot(X, E, linewidth = 3, color = 'k')
    plt.axhline(y=fermi-vacuum, linewidth = 2, color='r', linestyle='--')


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

        x[iz] = iz * (1 / ang2bohr * cell[2,2] / mesh[2])
        y[iz] = plane / count

    return x, y

files = sys.argv[1]
here = os.getcwd

cell, mesh, grid = io.readGrid(files)

x, y = planeaverage(cell, grid)

set_plot()   
plot_PAV(x/ang2bohr,y * Ry2eV-vacuum)

plt.tight_layout()
plt.savefig('PAV.png')
plt.close()
            
