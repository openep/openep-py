OpenEP-Py
=========

.. start-description

.. start-badges

|docs|
|testing|
|codecov|

.. |docs| image:: https://readthedocs.org/projects/openep-py/badge/?style=flat
    :target: https://openep-py.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |testing| image:: https://github.com/openep/openep-py/actions/workflows/python-app.yml/badge.svg
    :target: https://github.com/openep/openep-py/actions
    :alt: Testing Status

.. |codecov| image:: https://codecov.io/gh/openep/openep-py/branch/dev/graph/badge.svg
    :target: https://codecov.io/gh/openep/openep-py
    :alt: Coverage Status

.. end-badges

The open source solution for electrophysiology data analysis with Python.

Warning!
========

This software is in the early stages of development. As such, it is likely to have breaking changes introduced frequently. We will remove this warning once there is a stable API for the software.

Installation
============
Using git
---------

Install this library from a git clone: ::

    git clone https://github.com/openep/openep-py.git
    cd openep-py
    python3 -m pip install -e .


By installing OpenEP with the above command, you can update to the latest version anytime by going to the `openep-py` directory and typing `git pull`. However, if the the list of requirements changes, you will then need to type :code:`python3 -m pip install -e .` again in order to install the latest requirements.

Using PyPI
----------
We have not yet released a version of OpenEP on PyPI, but we're working on it.

.. end-description

Documentation
=============

The full documentation of OpenEP-Py is available at: `openep-py.readthedocs.io <https://openep-py.readthedocs.io/en/latest/?badge=latest>`__
