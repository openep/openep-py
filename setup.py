import os
from setuptools import setup, find_packages
from pkg_resources import parse_requirements

source_dir = os.path.abspath(os.path.dirname(__file__))

# read the version and other strings from _version.py
with open(os.path.join(source_dir, "openep/_version.py")) as o:
    exec(o.read())

# read install requirements from requirements.txt
with open(os.path.join(source_dir, "requirements.txt")) as o:
    requirements = [str(r) for r in parse_requirements(o.read())]

setup(
    name=__appname__,
    version=__version__,
    description=__description__,
    author=__author__,
    author_email=__author_email__,
    packages=find_packages(),
    install_requires=requirements,
    entry_points={},
)
