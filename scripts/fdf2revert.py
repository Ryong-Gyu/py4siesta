#!/usr/bin/env python
from NanoCore import *
import sys
import numpy as np

fname = sys.argv[1]
plane = sys.argv[2]

atom = s2.read_fdf(fname)
atom_new = atom.get_mirrored_structure(plane=plane)
sim = s2.Siesta(atom_new)
sim.write_struct()

