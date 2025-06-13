import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord
import astropy.units as u

from ..config import NSIDE

class Overlay():
    def __init__(self, coord: SkyCoord, color: str, radius: u.Quantity) -> "Overlay":
        if not radius.unit.is_equivalent(u.rad):
            raise ValueError("radius must have angular units")
        
        lat = coord.b.deg if coord.frame.name == "galactic" else coord.dec.deg
        long = coord.l.deg if coord.frame.name == "galactic" else coord.ra.deg


        theta = np.radians(90 - lat)
        phi = np.radians(long)
        obs_vec = hp.ang2vec(theta, phi)
        pixels = hp.query_disc(nside=NSIDE, vec=obs_vec, radius=np.radians(radius.to_value(u.deg)))
        theta, phi = hp.pix2ang(NSIDE, pixels)

        self.theta = theta
        self.phi = phi
        self.radius = radius.to_value(u.deg)
        self.color = color

    def _place_overlay(self):
        hp.projscatter(self.theta, self.phi, c=self.color, s=1)
