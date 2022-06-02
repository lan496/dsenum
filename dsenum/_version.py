# https://github.com/pypa/setuptools_scm/#retrieving-package-version-at-runtime
try:
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("dsenum")
except PackageNotFoundError:
    # package is not installed
    pass
