#!/bin/sh -eu
find -name "*.pyc" -delete
find -name "*.pyo" -delete
find -name "__pycache__" -delete
rm -r build dist
