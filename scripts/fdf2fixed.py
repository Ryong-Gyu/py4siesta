#!/usr/bin/env python
from NanoCore import *
import sys,os
import numpy as np

fname = sys.argv[1]
atom = s2.read_fdf(fname)
sim = s2.Siesta(atom)
sim.write_struct()
os.system(f'mv STRUCT.fdf {fname}')
