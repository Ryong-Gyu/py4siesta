#!/usr/bin/env python
from NanoCore import *
import sys
import numpy as np

#atom = s2.read_fdf('gnr.fdf')
atom1 = s2.read_fdf(sys.argv[1])
atom2 = s2.read_fdf(sys.argv[2])

atom3 = atom1.copy()

for initial_atom, final_atom, middle_atom in zip(atom1._atoms, atom2._atoms, atom3._atoms):

    initial_position = np.array(initial_atom.get_position(), dtype=float, copy=True)
    final_position = np.array(final_atom.get_position(), dtype=float, copy=True)
    middle_position = (initial_position + final_position) / 2
    middle_atom.set_position(Vector(middle_position))

sim = s2.Siesta(atom3)
sim.write_struct()

