# OpenEP-Py

[![Documentation Status](https://readthedocs.org/projects/openep-py/badge/?version=latest)](https://openep-py.readthedocs.io/en/latest/?badge=latest)
[![Testing](https://github.com/openep/openep-py/actions/workflows/python-app.yml/badge.svg)](https://github.com/openep/openep-py/actions)
[![codecov](https://codecov.io/gh/openep/openep-py/branch/dev/graph/badge.svg)](https://codecov.io/gh/openep/openep-py)

The open source solution for electrophysiology data analysis with Python.

# Warning!
This software is in the early stages of development. As such, it is likely to have breaking changes introduced frequently. We will remove this warning once there is a stable API for the software.

# Installation
## Using git
Install this library from a git clone:

```bash
git clone https://github.com/openep/openep-py.git
cd openep-py
python3 -m pip install -e .
```

By installing OpenEP with the above command, you can update to the latest version anytime by going to the `openep-py` directory and typing `git pull`. However, if the the list of requirements changes, you will then need to type `python3 -m pip install -e .` again in order to install the latest requirements.

## Using PyPI
We have not yet relased a version of OpenEP on PyPI, but we're working on it.
