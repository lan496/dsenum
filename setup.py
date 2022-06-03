from glob import glob
from os.path import basename, splitext
from pathlib import Path

from setuptools import Extension, find_packages, setup
from Cython.Build import cythonize

ext_modules = [Extension("dsenum.core", ["dsenum/core.pyx"], extra_compile_args=["-O3"])]

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
    python_requires=">=3.8",
    setup_requires=["setuptools_scm", "numpy", "Cython>=0.29.1"],
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
    extras_require={
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
            "sphinx-autobuild",
            "furo",
            "m2r2",
        ],
    },
    tests_require=["pytest"],
    ext_modules=cythonize(ext_modules),
    include_package_data=True,
    zip_safe=False,
    use_scm_version=True,
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
