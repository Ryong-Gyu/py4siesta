#!/usr/bin/env python
from NanoCore import *
import sys

fname = sys.argv[1]
direction = sys.argv[2]

name = fname.split('.')[0]


fdf = s2.read_fdf(fname)
fdf.select_all()
fdf.sort(option = direction)
fdf.set_serials(1)
sim = s2.Siesta(fdf)
sim.write_struct()
#io.write_xyz(name+'_new.xyz', xyz)
