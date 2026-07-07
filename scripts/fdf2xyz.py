from NanoCore import *
import sys

fname = sys.argv[1]

fdf = s2.read_fdf(fname)
io.write_xyz('STRUCT.xyz', fdf, comm=None, append=False)

