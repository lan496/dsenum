# dsenum

[![testing](https://github.com/lan496/dsenum/actions/workflows/testing.yml/badge.svg?branch=master)](https://github.com/lan496/dsenum/actions/workflows/testing.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/lan496/dsenum/master.svg)](https://results.pre-commit.ci/latest/github/lan496/dsenum/master)
[![codecov](https://codecov.io/gh/lan496/dsenum/branch/master/graph/badge.svg)](https://codecov.io/gh/lan496/dsenum)
[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/dsenum)
[![PyPI version](https://badge.fury.io/py/dsenum.svg)](https://badge.fury.io/py/dsenum)
![PyPI - Downloads](https://img.shields.io/pypi/dm/dsenum)

Derivative structure enumerator for multilattice

- Github: https://github.com/lan496/dsenum
- PyPI: https://pypi.org/project/dsenum/

## Installation

- Python>=3.8

Install in local:
```shell
git clone git@github.com:lan496/dsenum.git
cd dsenum
pip install -e .
```

For development:
```shell
git clone git@github.com:lan496/dsenum.git
cd dsenum
pip install -e ".[dev,docs]"
pre-commit install
```

## Usage

```python
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie

from dsenum import StructureEnumerator

latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
coords = [[0, 0, 0]]
structure = Structure(latt, [DummySpecie('X')] * len(coords), coords)

num_type = 2
index = 4

dstructs = StructureEnumerator(structure, index, num_type).generate()
print(len(dstructs))  # -> 12
```

See `examples/Sn_oxide.py` for more complicated use cases.

## Official Implementation
- https://github.com/msg-byu/enumlib

## References
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.224115
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.80.014120
- https://iopscience.iop.org/article/10.1088/0953-8984/25/10/105401
- http://blog.dlfer.xyz/post/2016-10-27-smith-normal-form/
