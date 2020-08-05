from os.path import splitext, basename, abspath, dirname, join
from setuptools import setup, Extension, find_packages
from Cython.Distutils import build_ext
from glob import glob

ext_modules = [Extension("dsenum.core", sources=["dsenum/core.pyx"], extra_compile_args=["-O3"])]


def get_version(rel_path):
    here = abspath(dirname(__file__))
    with open(join(here, rel_path), "r") as f:
        for line in f.read().splitlines():
            if line.startswith("__version__"):
                # e.g. __version__ = "0.1.0"
                delim = '"' if line else "'"
                return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")


setup(
    name="dsenum",
    version=get_version("dsenum/__init__.py"),
    license="MIT",
    description="Derivative structure enumerator for multilattice",
    # long_description="",
    author="Kohei Shinohara",
    author_email="kohei19950508@gmail.com",
    packages=find_packages("dsneum"),
    py_modules=[splitext(basename(path))[0] for path in glob("dsenum/*.py")],
    python_requires=">=3.6",
    install_requires=["setuptools", "Cython"],
    tests_require=["pytest"],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    include_package_data=True,
    # extras_requires={},
)
