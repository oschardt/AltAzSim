import astropy.units as u
from astropy.units import Quantity
from astropy.coordinates import EarthLocation

import numpy as np

class Telescope:
    def __init__(self, location: EarthLocation, az_slew: Quantity, tip_slew: Quantity, az_accel: Quantity, tip_accel: Quantity,
                 max_alt: Quantity, min_alt: Quantity, max_az: Quantity, min_az: Quantity, angular_res: Quantity):
        # make sure variables are of the right type
        if not isinstance(location, EarthLocation):
            raise ValueError("location must be an astropy.coordinates.EarthLocation instance")
        if not az_slew.unit.is_equivalent(u.deg / u.s):
            raise ValueError("az_slew must have astropy angular velocity units")
        if not None and not az_accel.unit.is_equivalent(u.deg / u.s ** 2):
            raise ValueError("az_accel must have astropy angular acceleration units")
        if not tip_slew.unit.is_equivalent(u.deg / u.s):
            raise ValueError("tip_slew must have astropy angular velocity units")
        if not None and not tip_accel.unit.is_equivalent(u.deg / u.s ** 2):
            raise ValueError("tip_accel must have astropy angular acceleration units")
        if not max_alt.unit.is_equivalent(u.deg):
            raise ValueError("max_alt must have astropy angular units")
        if not min_alt.unit.is_equivalent(u.deg):
            raise ValueError("min_alt must have astropy angular units")
        if not max_az.unit.is_equivalent(u.deg):
            raise ValueError("max_az must have astropy angular units")
        if not min_az.unit.is_equivalent(u.deg):
            raise ValueError("min_az must have astropy angular units")
        if not angular_res.unit.is_equivalent(u.deg**2):
            raise ValueError("angular_res_az must have astropy solid angle units")
        
        # initialize variables
        self.location = location
        
        # everything internally uses just the value
        self._az_slew = az_slew.to_value(u.deg / u.s)
        self._az_accel = az_accel.to_value(u.deg/u.s**2)
        self._tip_slew = tip_slew.to_value(u.deg / u.s)
        self._tip_accel = tip_accel.to_value(u.deg/u.s**2)
        self._max_alt = max_alt.to_value(u.deg)
        self._min_alt = min_alt.to_value(u.deg)
        self._max_az = max_az.to_value(u.deg)
        self._min_az = min_az.to_value(u.deg)
        self._angular_res = angular_res.to_value(u.deg**2)

        self._hp_query_radius = self._get_query_radius()
    
    def _get_query_radius(self):
        omega = self._angular_res * (np.pi / 180)**2
        return np.arccos(1 - omega / (2 * np.pi))

    @property
    def az_slew(self):
        return self._az_slew * (u.deg / u.s)

    @az_slew.setter
    def az_slew(self, val):
        if not val.unit.is_equivalent(u.deg / u.s):
            raise ValueError("az_slew must have astropy angular velocity units")
        self._az_slew = val.to_value(u.deg / u.s)

    @property
    def az_accel(self):
        if self._az_accel is None:
            raise ValueError("az_accel is not initialized")
        return self._az_accel * (u.deg / u.s**2)

    @az_accel.setter
    def az_accel(self, val):
        if not val.unit.is_equivalent(u.deg / u.s**2):
            raise ValueError("az_accel must have astropy angular acceleration units")
        self._az_accel = val.to_value(u.deg / u.s**2)

    @property
    def tip_slew(self):
        return self._tip_slew * (u.deg / u.s)

    @tip_slew.setter
    def tip_slew(self, val):
        if not val.unit.is_equivalent(u.deg / u.s):
            raise ValueError("tip_slew must have astropy angular velocity units")
        self._tip_slew = val.to_value(u.deg / u.s)

    @property
    def tip_accel(self):
        if self._tip_accel is None:
            raise ValueError("tip_accel is not initialized")
        return self._tip_accel * (u.deg / u.s**2)

    @tip_accel.setter
    def tip_accel(self, val):
        if not val.unit.is_equivalent(u.deg / u.s**2):
            raise ValueError("tip_accel must have astropy angular acceleration units")
        self._tip_accel = val.to_value(u.deg / u.s**2)

    @property
    def max_alt(self):
        return self._max_alt * u.deg

    @max_alt.setter
    def max_alt(self, val):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("max_alt must have astropy angular units")
        self._max_alt = val.to_value(u.deg)

    @property
    def min_alt(self):
        return self._min_alt * u.deg

    @min_alt.setter
    def min_alt(self, val):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("min_alt must have astropy angular units")
        self._min_alt = val.to_value(u.deg)

    @property
    def max_az(self):
        return self._max_az * u.deg

    @max_az.setter
    def max_az(self, val):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("max_az must have astropy angular units")
        self._max_az = val.to_value(u.deg)

    @property
    def min_az(self):
        return self._min_az * u.deg

    @min_az.setter
    def min_az(self, val):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("min_az must have astropy angular units")
        self._min_az = val.to_value(u.deg)

    @property
    def angular_res(self):
        return self._angular_res * u.deg**2
    
    @angular_res.setter
    def angular_res(self, val):
        if not val.unit.is_equivalent(u.deg**2):
            raise ValueError("angular_res must have astropy solid angle units")
        self._min_az = val.to_value(u.deg**2)


