from setuptools import setup, Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension('dsenum.core', sources=['dsenum/core.pyx'])
]

setup(
    packages=['dsenum'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest', 'pytest-cov'],
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext}
)
