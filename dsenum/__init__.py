import os
import subprocess

from ._version import get_versions  # type: ignore

from dsenum.enumerate import (
    StructureEnumerator,
    ZddStructureEnumerator,
)


__version__ = get_versions()["version"]
del get_versions


def get_git_commit_hash() -> str:
    cwd = os.path.dirname(os.path.abspath(__file__))
    out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=cwd)
    return out.strip().decode("ascii")


def get_version() -> str:
    out = __version__ + "+" + get_git_commit_hash()
    return out
