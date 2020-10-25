# pyzdd
Python wrapper to TdZdd

## Installation

```script
git clone --recursive git@github.com:lan496/pyzdd.git
```

### Conda
```
conda create --name pyzdd python=3.8 pip
conda activate pyzdd
pip install -r requirements.txt
pip install -e .
```

## Development

### Installation
```
conda create --name pyzdd python=3.8 pip
pyenv activate pyzdd
./clean.sh
pip install -r requirements.txt
pip install -r requirements-dev.txt
pip install -e .
pre-commit install
```

### Write Custom Specification
1. Write a TdZdd-specification in `src/spec/*.hpp`
2. Let the new specification class be `A`, wrap the following classes and methods in `src/wrapper.cpp`
    - `tdzdd::DdSpecBase<A, 2>`
    - `tdzdd::DdSpec<A, T, 2>`
    - `A`
    - `tdzdd::DdStructure<2>::zddSubset<A>`
3. import `_pyzdd.A` in `pyzdd/__init__.py`

## References
- https://github.com/kunisura/TdZdd
