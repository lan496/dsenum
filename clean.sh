#!/bin/bash -u
find -name "*.so" -delete
find -name "*.pyc" -delete
find -name "*.pyo" -delete
find -name "__pycache__" -delete
rm -rf build dist *.egg-info
