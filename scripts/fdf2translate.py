#!/usr/bin/env python
from NanoCore import *
import sys
import numpy as np

fname = sys.argv[1]
x = float(sys.argv[2])
y = float(sys.argv[3])
z = float(sys.argv[4])

atom = s2.read_fdf(fname)
atom.select_all()
distance = np.array([x,y,z])
atom.translate(*distance)

sim = s2.Siesta(atom)
sim.write_struct()
