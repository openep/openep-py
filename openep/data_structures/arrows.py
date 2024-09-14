from attr import attrs
import numpy as np

__all__ = []


@attrs(auto_attribs=True, auto_detect=True)
class Arrows:
    """
    Class for storing information about arrows and lines on surface

    Args:
        fibres (np.ndarray): array of shape N_cells x 3
        divergence (np.ndarray): array of shape N_cells x 3
        linear_connections (np.ndarray): array of shape M x 3 (represents the linear connections between endo and epi)
        linear_connection_regions (np.ndarray): array of shape N_cells
    """

    # TODO: move divergence arrows into Arrows class
    # TODO: remove longitudinal and transversal arrows from Fields class
    fibres: np.ndarray = None
    divergence: np.ndarray = None
    linear_connections: np.ndarray = None
    linear_connection_regions: np.ndarray = None

    def __repr__(self):
        return f"arrows: {tuple(self.__dict__.keys())}"

    def __getitem__(self, arrow):
        try:
            return self.__dict__[arrow]
        except KeyError:
            raise ValueError(f"There is no arrow '{arrow}'.")

    def __setitem__(self, arrow, value):
        if arrow not in self.__dict__.keys():
            raise ValueError(f"'{arrow}' is not a valid arrow name.")
        self.__dict__[arrow] = value

    def __iter__(self):
        return iter(self.__dict__.keys())

    def __contains__(self, arrow):
        return arrow in self.__dict__.keys()

    @property
    def linear_connection_regions_names(self):
        if self.linear_connection_regions is None:
            return []
        regions = np.unique(self.linear_connection_regions).astype(str)
        return regions.tolist()

    def copy(self):
        """Create a deep copy of Arrows"""

        arrows = Arrows()
        for arrow in self:
            if self[arrow] is None:
                continue
            arrows[arrow] = np.array(self[arrow])

        return arrows
