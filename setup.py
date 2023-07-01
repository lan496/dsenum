from pathlib import Path

from setuptools import Extension, setup
from Cython.Build import cythonize

ext_modules = [Extension("dsenum.core", ["src/dsenum/core.pyx"], extra_compile_args=["-O3"])]

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
    package_dir={"": "src"},
    package_data={"dsenum": ["py.typed"]},
    python_requires=">=3.8",
    setup_requires=["setuptools_scm", "numpy", "Cython>=0.29.1"],
    install_requires=[
        "setuptools",
        "setuptools_scm",
        "wheel",
        "pymatgen>=2021.3.9",
        # Blocked until next version of pymatgen>2022.11.7 is released
        "numpy<1.26",
        "sympy",  # sympy.utilities.iterables.multiset_permutations
        "scipy",  # scipy.special.binom
        "tqdm",
        "hsnf>=0.3.15",
        "pyzdd==0.2.6",
    ],
    extras_require={
        "dev": [
            "ipython",
            "notebook",
            "jupyter_contrib_nbextensions",
            "matplotlib",
            "seaborn",
            "pandas",
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
            "nbsphinx==0.9.2",
            "myst-parser",
            "sphinx-book-theme",
            "ipykernel",
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
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
