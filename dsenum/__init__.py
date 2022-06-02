from dsenum.enumerate import StructureEnumerator

# https://github.com/pypa/setuptools_scm/#retrieving-package-version-at-runtime
try:
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("hsnf")
except PackageNotFoundError:
    # package is not installed
    pass
