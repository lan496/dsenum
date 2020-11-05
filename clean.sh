#!/bin/sh -eu
BASEDIR=$(dirname "$0")

find ${BASEDIR} -name "*.pyc" -delete
find ${BASEDIR} -name "*.pyo" -delete
find ${BASEDIR} -name "__pycache__" -delete
rm -rf "${BASEDIR}/build" "${BASEDIR}/dist"
