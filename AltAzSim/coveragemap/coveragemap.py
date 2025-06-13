import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from ..config import NSIDE

from .overlay import Overlay


class CoverageMap():
    def __init__(self) -> "CoverageMap":
        self.nside = NSIDE
        self.map = self._new_map(NSIDE)
        self.title = "Scans Coverage"
        self.unit = "Passes"

        # used for displaying
        self.norm = "norm"
        self.show_title = True
        self.show_unit_bar = True
        self.threshold = None

    @classmethod
    def from_map(cls, map: np.ndarray, title: str = "Pre-existing Map") -> "CoverageMap":
        nside = hp.get_nside(map)
        new_map = cls()
        new_map.nside = nside
        new_map.title = title
        new_map.map = map
        return new_map


    def copy(self):
        copy = CoverageMap()
        copy.nside = self.nside
        copy.map = self.map.copy()
        copy.title = self.title
        copy.unit = self.unit
        return copy

    def _new_map(self, nside):
        if not(nside > 0 and (nside & (nside - 1)) == 0):
            raise ValueError("nside must be a power of 2")
        npix = hp.nside2npix(nside)
        return np.zeros(npix)
    
    def _prep_for_disp(self):
        """
        hp.mollview()'s based on preferences
        """

        title = self.title if self.show_title else None
        map = self.map.copy()
        map[map < 1] = hp.UNSEEN
        
        if self.threshold is not None:
            mask = map >= self.threshold
            map[mask] = self.threshold
        else:
            # norm = "norm"
            alpha = None

        if self.show_unit_bar:
            cbar = True
            unit = self.unit
        else:
            cbar = False
            unit = None
        
        hp.mollview(map, title=title, unit=unit, cbar=cbar, cmap=plt.cm.viridis, norm=self.norm)
        
    def set_display_pref(self, show_title: bool = True, show_unit_bar: bool = True, threshold: int = None, norm: str = "norm"):
        if not isinstance(show_title, bool):
            raise ValueError("show_title must be bool")
        if not isinstance(show_unit_bar, bool):
            raise ValueError("show_unit_bar must be bool")
        if threshold is not None and not isinstance(threshold, int):
            raise ValueError("threshold must be int or None")
        
        self.show_title = show_title
        self.show_unit_bar = show_unit_bar
        self.threshold = threshold
        self.norm = norm

    def display(self, overlays: None | Overlay | list[Overlay] = None):
        """
        Displays the map

        Parameters
        -------
        overlays  :  Overlay or list of Overlay, optional
          points to be overlayed on the map
        """
        if (overlays is not None) and (not isinstance(overlays, list)):
            overlays = [overlays]

        self._prep_for_disp()

        if overlays is not None:
            for overlay in overlays:
                overlay._place_overlay()

        plt.show()
        plt.close()

    def get_coverage(self) -> np.floating:
        '''
        Returns the fraction of pixels that were covered
        '''
        return np.mean(self.map >= 1)
