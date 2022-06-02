# https://github.com/pypa/setuptools_scm/#retrieving-package-version-at-runtime
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("dsenum")
except PackageNotFoundError:
    # package is not installed
    pass
