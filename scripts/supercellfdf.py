#!/usr/bin/env python
from NanoCore import *
import sys

fname = sys.argv[1]
x = sys.argv[2]
y = sys.argv[3]
z = sys.argv[4]

atom = s2.read_fdf(fname)
atom2 = atom * [int(x), int(y), int(z)]
sys = s2.Siesta(atom2)
sys.write_struct()
