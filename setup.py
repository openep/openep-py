import os
import glob
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
    py_modules=[os.path.splitext(os.path.basename(path))[0] for path in glob.glob('openep/*.py')],
    install_requires=requirements,
    url='https://github.com/openep/openep-gui',
    project_urls={
        'Documentation': 'https://openep-py.readthedocs.io/en/latest/',
        'Issue Tracker': 'https://github.com/openep/openep-py/issues',
    },
    python_requires='>=3.8',
    include_package_data=True,
    package_data={
        'openep/_datasets/OpenEP-MATLAB': ['_datasets/OpenEP-MATLAB/*.mat']
    }
)
