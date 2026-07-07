#!/usr/bin/env python
from NanoCore import *
import sys, os
import glob

files = glob.glob('*STRUCT_OUT')[0]

atom = s2.read_struct_out(files)
system = s2.Siesta(atom)
system.write_struct()
