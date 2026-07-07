#!/usr/bin/env python
from NanoCore import *
import sys

fdf1 = sys.argv[1]
fdf2 = sys.argv[2]

atm1 = s2.read_fdf(fname)
atm2 = s2.read_fdf(fname)

atom3 = atm1 + atm2

sim = s2.Siesta(atom3)
sim.write_struct()
#io.write_xyz(name+'_new.xyz', xyz)
