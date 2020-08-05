# dsenum
![Build Status](https://travis-ci.com/lan496/dsenum.svg?branch=master)
[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)
[![Maintainability](https://api.codeclimate.com/v1/badges/ba967883e004590a2a23/maintainability)](https://codeclimate.com/github/lan496/dsenum/maintainability)
[![Test Coverage](https://api.codeclimate.com/v1/badges/ba967883e004590a2a23/test_coverage)](https://codeclimate.com/github/lan496/dsenum/test_coverage)

Derivative structure enumerator for multilattice


## installation

### Requirements
- cython

### install

#### pip
```
git clone git@github.com:lan496/dsenum.git
cd dsenum
pip install Cython
pip install -r requirements.txt
python setup.py install
pre-commit install
```

#### conda
```script
git clone git@github.com:lan496/dsenum.git
cd dsenum
conda env create -f environment.yml
conda activate dsenum
pip install -e .
pre-commit install
```

## Usage

```sample.py
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie

from dsenum.enumerate import enumerate_derivative_structures

latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
coords = [[0, 0, 0]]
structure = Structure(latt, [DummySpecie('X')] * len(coords), coords)

num_type = 2
index = 4

dstructs = enumerate_derivative_structures(structure, index, num_type)
print(len(dstructs))  # -> 12
```

## Development

### Github Actions
To use `actions/setup-python`, we need to use the full container (costs >18GB...)

```
act -P ubuntu-latest=nektos/act-environments-ubuntu:18.04
```

### Pre-commit
- black
    see `pyproject.toml` and `.pre-commit-config.yaml`
    if you do not want to format codes, wrap the block by `fmt: off/on`
    ```
    # fmt: off
    something_you_do_not_want_to_be_formatted
    # fmt: on
    ```

## Official Implementation
- https://github.com/msg-byu/enumlib

## References
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.224115
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.80.014120
- https://iopscience.iop.org/article/10.1088/0953-8984/25/10/105401
- http://blog.dlfer.xyz/post/2016-10-27-smith-normal-form/
