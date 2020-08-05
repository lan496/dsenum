#!/bin/bash -eu
find -name "*.pyc" -delete
find -name "*.pyo" -delete
find -name "__pycache__" -delete
rm -r build
