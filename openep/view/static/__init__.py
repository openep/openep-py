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

# Icons for the link/unlink views button of secondary viewers
LINK_ICON = resource_filename(__name__,
                              "buttons/views-linked.png")
UNLINK_ICON = resource_filename(__name__,
                                "buttons/views-unlinked.png")
