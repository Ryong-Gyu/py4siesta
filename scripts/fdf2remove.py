#!/usr/bin/env python
from NanoCore import *
import sys
import numpy as np

atom = s2.read_fdf('STRUCT.fdf')
atom.select_elements('C')
atom.delete()
sim = s2.Siesta(atom)
sim.write_struct()
