#!/usr/bin/env python
from NanoCore import *
import sys
    
fdf_name = sys.argv[1]
at = s2.read_fdf(fdf_name)
vasp.write_poscar(at, 'POSCAR')
