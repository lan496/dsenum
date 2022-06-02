from glob import glob
from os.path import basename, splitext
from pathlib import Path

from setuptools import Extension, find_packages, setup

ext_modules = [Extension("dsenum.core", sources=["dsenum/core.c"], extra_compile_args=["-O3"])]

# Import the README and use it as the long-description.
with open(Path(__file__).resolve().parent / "README.md") as f:
    long_description = "\n" + f.read()


setup(
    name="dsenum",
    license="MIT",
    description="Derivative structure enumerator for multilattice",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Kohei Shinohara",
    author_email="kohei19950508@gmail.com",
    packages=find_packages("dsneum"),
    package_data={"dsenum": ["py.typed"]},
    py_modules=[splitext(basename(path))[0] for path in glob("dsenum/*.py")],
    python_requires=">=3.7",
    setup_requires=["setuptools_scm", "numpy"],
    install_requires=[
        "setuptools",
        "setuptools_scm",
        "wheel",
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
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
