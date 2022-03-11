import os
from setuptools import setup, find_packages
from pkg_resources import parse_requirements

source_dir = os.path.abspath(os.path.dirname(__file__))

# read the version and other strings from _version.py
version_info = {}
with open(os.path.join(source_dir, "openep/_version.py")) as o:
    exec(o.read(), version_info)

# read install requirements from requirements.txt
with open(os.path.join(source_dir, "requirements.txt")) as o:
    requirements = [str(r) for r in parse_requirements(o.read())]

setup(
    name="OpenEp",
    version=version_info['__version__'],
    description="Open Source solution for electrophysiology data analysis",
    author="Steven Williams",
    author_email="steven.williams@ed.ac.uk",
    packages=find_packages(),
    install_requires=requirements,
)
