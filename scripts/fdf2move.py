#!/usr/bin/env python
from NanoCore import *
import sys
import numpy as np

#atom = s2.read_fdf('gnr.fdf')
atom = s2.read_fdf('base.fdf')
#atom = s2.read_fdf('backbone.fdf')

atom.select_all()
atom.translate(0,0.25,0)

'''
atom.select_atmnbs(list(range(1,369)))
atom.translate(0,-0.799200000,0)

atom.select_atmnbs(list(range(369,721)))
atom.translate(0,0.7744,0)
'''

sim = s2.Siesta(atom)
sim.write_struct()

