#!/usr/bin/env python
from pymatgen.core import Structure
import sys

structure = Structure.from_file(sys.argv[1])
primitive_structure = structure.get_primitive_structure()
primitive_structure.to(fmt="poscar", filename="POSCAR_primitive")
