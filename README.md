# dsenum
![Build Status](https://travis-ci.com/lan496/dsenum.svg?branch=master)
[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

Derivative structure enumerator for multilattice

## install
```
git clone git@github.com:lan496/dsenum.git
cd dsenum
pip install -r requirements.txt
pip install -e .
```

## Usage

```sample.py
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie

from dsenum.enumerate import enumerate_derivative_structures

latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
coords = [[0, 0, 0]]
struct = Structure(latt, [DummySpecie('X')] * len(coords), coords)

num_type = 2
index = 4

dstructs = enumerate_derivative_structures(structure, index, num_type)
print(len(dstructs))  # ->14
```

## Official Implementation
- https://github.com/msg-byu/enumlib

## References
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.224115
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.80.014120
- https://iopscience.iop.org/article/10.1088/0953-8984/25/10/105401
- http://blog.dlfer.xyz/post/2016-10-27-smith-normal-form/
