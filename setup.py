from glob import glob
from os.path import basename, splitext

from setuptools import Extension, find_packages, setup

import versioneer

ext_modules = [Extension("dsenum.core", sources=["dsenum/core.c"], extra_compile_args=["-O3"])]


setup(
    name="dsenum",
    version=versioneer.get_version(),  # type: ignore
    license="MIT",
    description="Derivative structure enumerator for multilattice",
    # long_description="",
    author="Kohei Shinohara",
    author_email="kohei19950508@gmail.com",
    packages=find_packages("dsneum"),
    py_modules=[splitext(basename(path))[0] for path in glob("dsenum/*.py")],
    python_requires=">=3.7",
    install_requires=[
        "setuptools",
        "pymatgen>=2020.4.29",
        "numpy",
        "sympy",
        "scipy",
        "matplotlib",
        "joblib",
        "tqdm",
    ],
    extras_requires={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov",
            "pytest-benchmark",
            "pre-commit",
            "black",
            "flake8",
            "mypy",
            "isort",
            "pyupgrade",
            "versioneer",
        ],
        "docs": [
            "sphinx",
            "sphinx-rtd-theme",
        ],
    },
    tests_require=["pytest"],
    ext_modules=ext_modules,
    include_package_data=True,
    zip_safe=False,
)
