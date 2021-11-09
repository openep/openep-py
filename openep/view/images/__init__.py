"""
OpenEP Images for the GUI
=========================

"""
__all__ = [
    "LOGO",
]

from pkg_resources import resource_filename

LOGO = resource_filename(__name__,
                         "openep-logo.png")