"""The setup script."""
from typing import Dict

from setuptools import setup

version: Dict[str, str] = {}
with open("mwrpy/version.py") as fp:
    exec(fp.read(), version)

setup(version=version["version"])
