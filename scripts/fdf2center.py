#!/usr/bin/env python
from NanoCore import *
import sys
import numpy as np

fname = sys.argv[1]
atom = s2.read_fdf(fname)
atom.select_all()
center = atom.center(mode="geom")
cell = np.array(atom.get_cell())


distance = (cell[0]+cell[1]+cell[2])/2 - np.array(center)

atom.translate(*distance)

sim = s2.Siesta(atom)
sim.write_struct()
