#!/usr/bin/env python
from NanoCore import *
import sys
import numpy as np

fname1 = sys.argv[1]
fname2 = sys.argv[2]

atom1 = s2.read_fdf(fname1)
atom2 = s2.read_fdf(fname2)
atom = atom1 + atom2
sim = s2.Siesta(atom)
sim.write_struct()


