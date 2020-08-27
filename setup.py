from os.path import splitext, basename
from setuptools import setup, Extension, find_packages
from glob import glob

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
    install_requires=["setuptools"],
    tests_require=["pytest"],
    ext_modules=ext_modules,
    include_package_data=True,
    # extras_requires={},
    zip_safe=False,
)
