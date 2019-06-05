from setuptools import setup, Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension('dsenum.core', sources=['dsenum/core.pyx'])
]

setup(
    packages=['dsenum'],
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext}
)
