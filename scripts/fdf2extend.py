#!/usr/bin/env python
from NanoCore import *
import sys
import numpy as np

fname = sys.argv[1]
ratio = float(sys.argv[2])
axis = sys.argv[3]

atom = s2.read_fdf(fname)

if axis == 'a':
    atom2 = atom.adjust_cell_size(ratio, direction=1)
elif axis == 'b':
    atom2 = atom.adjust_cell_size(ratio, direction=2)
elif axis == 'c':
    atom2 = atom.adjust_cell_size(ratio, direction=3)

sim = s2.Siesta(atom2)
sim.write_struct()
