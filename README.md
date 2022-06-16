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

dsenum works with Python3.8+ and can be installed via PyPI:

```shell
pip install dsenum
```

Or in local:
```shell
git clone git@github.com:lan496/dsenum.git
cd dsenum
pip install -e .
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

See [docs/examples/constraints.ipynb](docs/examples/constraints.ipynb) for more complicated use cases.

## How to cite dsenum

If you use `dsenum.ZddStructureEnumerator` in your research, please cite the following articles.

```
@article{doi:10.1063/5.0021663,
    author = {Shinohara,Kohei  and Seko,Atsuto  and Horiyama,Takashi  and Ishihata,Masakazu  and Honda,Junya  and Tanaka,Isao },
    title = {Enumeration of nonequivalent substitutional structures using advanced data structure of binary decision diagram},
    journal = {J. Chem. Phys.},
    volume = {153},
    number = {10},
    pages = {104109},
    year = {2020},
    doi = {10.1063/5.0021663},
    URL = {https://doi.org/10.1063/5.0021663},
}
```

```
@inproceedings{Horiyama2018,
  memo ={Isomorphism Elimination by Zero-Suppressed Binary Decision Diagrams},
  author={Takashi Horiyama and Masahiro Miyasaka and Riku Sasaki},
  booktitle={the Canadian Conference on Computational Geometry},
  pages={360--366},
  address={Winnipeg, Manitoba, Canada}
  year={2018},
  url={http://www.cs.umanitoba.ca/~cccg2018/papers/session7B-p2.pdf}
}
```

## Acknowledgments

I acknowledge Dr. Takashi Horiyama for sharing his implementation of the frontier method for isomorphism-elimination decision diagram.
I also appreciate his kindness to allow publishing the code.
