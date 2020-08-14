# setup.py reflects the below version string.
__version__ = "0.3.2"

import os
import subprocess

from dsenum.enumerate import StructureEnumerator


def get_git_commit_hash() -> str:
    cwd = os.path.dirname(os.path.abspath(__file__))
    out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=cwd)
    return out.strip().decode("ascii")


def get_version() -> str:
    out = __version__ + "+" + get_git_commit_hash()
    return out
