## How to compile and check documents

```
cd docs
sphinx-apidoc -f -o source ../dsenum
make html
```

```
cd build/html
python -m http.server
```
