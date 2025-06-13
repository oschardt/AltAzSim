from ..utils.linearmotion import linear_motion, linear_motion_dur, xf_linear_motion_from_time
from ..utils import format_seconds
from .scanprofile import ScanProfile
from ..telescope import Telescope

import astropy.units as u
from astropy.units import Quantity
import numpy as np

class BackAndForthScan(ScanProfile):
    def __init__(self, *args):
        raise RuntimeError("Use for_duration() or for_count() to construct this class.")
    
    _DEFAULT_CONFIG = {
        "include_alignment_time": True,
        "show_alignment": True,
    }
    
    def _init_internal(self, scan_mode: str, max_angle: Quantity, min_angle: Quantity, slew: Quantity, fixed_deg: Quantity, count: int | None = None, duration: Quantity | None = None, config: dict | None = None):
        super()._init_internal(duration, count)
        
        if not isinstance(scan_mode, str):
            raise ValueError("scan_mode must be a string")
        if scan_mode not in ["az", "tip"]:
            raise ValueError("scan_mode must be either 'az' or 'tip'")
        
        # make sure variables are of the right type
        if not max_angle.unit.is_equivalent(u.deg):
            raise ValueError("max_angle must have angular units")
        if not min_angle.unit.is_equivalent(u.deg):
            raise ValueError("min_angle must have angular units")
        if not slew.unit.is_equivalent(u.deg/u.s):
            raise ValueError("slew must have angular velocity units")
        if not fixed_deg.unit.is_equivalent(u.deg):
            raise ValueError("fixed_deg must have angular units")
        if count is not None and not isinstance(count, int):
            raise ValueError("count must be int")
        if duration is not None and not duration.unit.is_equivalent(u.s):
            raise ValueError("duration must have time units")
        
        self._scan_mode = scan_mode
        self._max_angle = max_angle.to_value(u.deg)
        self._min_angle = min_angle.to_value(u.deg)
        self._slew = slew.to_value(u.deg/u.s)
        self._fixed_deg = fixed_deg.to_value(u.deg)

        if config is None:
            self._config = self._DEFAULT_CONFIG.copy()
        else:
            invalid_keys = set(config) - self._DEFAULT_CONFIG.keys()
            if invalid_keys:
                raise KeyError(f"Invalid config keys: {invalid_keys}")
            self._config = config

            

    
    @classmethod
    def from_duration(cls, scan_mode: str, max_angle: Quantity, min_angle: Quantity, slew: Quantity, fixed_deg: Quantity, duration: Quantity, config: dict = None):
        self = cls.__new__(cls)
        self._init_internal(scan_mode, max_angle, min_angle, slew, fixed_deg, None, duration, config)
        return self
    
    @classmethod
    def from_count(cls, scan_mode: str, max_angle: Quantity, min_angle: Quantity, slew: Quantity, fixed_deg: Quantity, count: Quantity, config: dict = None):
        self = cls.__new__(cls)
        self._init_internal(scan_mode, max_angle, min_angle, slew, fixed_deg, count, None, config)
        return self

    def copy(self):
        copy = self.__class__.__new__(self.__class__)
        copy._init_internal(
            self._scan_mode,
            self._max_angle * u.deg,
            self._min_angle * u.deg,
            self._slew * u.deg/u.s,
            self._fixed_deg * u.deg,
            self._count if not self._is_dur_based else None,
            self._duration * u.s if self._is_dur_based else None,
            self._config.copy()
        )
        return copy
    
    def scan_info(self):
        return (
            f"{"azimuth" if self._scan_mode == "az" else "tip"} range: {self._min_angle} degs to {self._max_angle} degs\n"
            f"slew: {self._slew:.4f} degs/s\n"
            f"fixed {"elevation" if self._scan_mode == "az" else "azimuth"}: {self._fixed_deg} degs\n"
        )
    
    @property
    def scan_mode(self) -> str:
        return self._scan_mode
    
    @scan_mode.setter
    def scan_mode(self, val):
        raise ValueError("scan_mode cannot be changed after initializaiton")


    @property
    def max_angle(self) -> Quantity:
        return self._max_angle * u.deg

    @max_angle.setter
    def max_angle(self, val: Quantity):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("max_angle must have angular units")
        self._max_angle = val.to_value(u.deg)


    @property
    def min_angle(self) -> Quantity:
        return self._min_angle * u.deg

    @min_angle.setter
    def min_angle(self, val: Quantity):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("min_angle must have angular units")
        self._min_angle = val.to_value(u.deg)


    @property
    def slew(self) -> Quantity:
        return self._slew * (u.deg / u.s)

    @slew.setter
    def slew(self, val: Quantity):
        if not val.unit.is_equivalent(u.deg / u.s):
            raise ValueError("slew must have angular velocity units")
        self._slew = val.to_value(u.deg / u.s)


    @property
    def fixed_deg(self) -> Quantity:
        return self._fixed_deg * u.deg

    @fixed_deg.setter
    def fixed_deg(self, val: Quantity):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("fixed_deg must have angular units")
        self._fixed_deg = val.to_value(u.deg)

    @classmethod
    def get_default_config(self) -> dict:
        return self._DEFAULT_CONFIG.copy()
    
    @property
    def config(self) -> dict:
        return self._config


    def validate_for_tel(self, telescope: Telescope):
        if self._scan_mode == "az" and (self._max_angle > telescope._max_az or self._min_angle < telescope._min_az):
            raise ValueError("Azimuth scanning angles out of range.")
        elif self._scan_mode == "tip" and (self._max_angle > telescope._max_alt or self._min_angle < telescope._min_alt):
            raise ValueError("Tip scanning angles out of range.")
        if self._scan_mode == "az" and (self._fixed_deg > telescope._max_alt or self._fixed_deg < telescope._min_alt):
            raise ValueError("Fixed degree of scan is out of azimuthal scanning range")
        elif self._scan_mode == "tip" and (self._fixed_deg > telescope._max_az or self._fixed_deg < telescope._min_az):
            raise ValueError("Fixed degree of scan is out of tip scanning range")
        if self._scan_mode == "az" and (self._slew > telescope._az_slew):
            raise ValueError("slew is not valid for telescope")
        if self._scan_mode == "tip" and (self._slew > telescope._tip_slew):
            raise ValueError("slew is not valid for telescope")


    def run(self):
        """
        Slews in one primary direction
        """
        is_az_scan = self._scan_mode == "az"
        prim_dir = self._init_az if is_az_scan else self._init_alt
        prim_slew = self._slew
        prim_accel = self._tel._az_accel if is_az_scan else self._tel._tip_accel
        fixed_dir = self._init_alt if is_az_scan else self._init_az
        fixed_slew = self._tel._tip_slew if is_az_scan else self._tel._az_slew
        fixed_accel = self._tel._tip_accel if is_az_scan else self._tel._az_accel

        time_elapsed = np.float64(0)
        
        prim_scan_vals = []
        fixed_scan_vals = []
        times = []

        # align secondary axis
        if fixed_dir != self._fixed_deg:
            if self._config["show_alignment"]:
                x,t = linear_motion(fixed_dir, self._fixed_deg, fixed_slew, fixed_accel, self._dt)
                y = np.full(len(x), prim_dir)

                fixed_scan_vals.append(x)
                prim_scan_vals.append(y)
                times.append(t)
            else:
                fixed_scan_vals.append([self._fixed_deg])
                prim_scan_vals.append([prim_dir])
                t = linear_motion_dur(fixed_dir, self._fixed_deg, fixed_slew, fixed_accel)
                times.append([t])

            if self._config["include_alignment_time"]:
                time_elapsed += times[-1][-1]
                if self._is_dur_based and time_elapsed > self._duration:
                    raise ValueError("scan does not allot enough time to begin scan. Secondary axis cannot align")

            fixed_dir = self._fixed_deg
        
        in_limits = (prim_dir > self._min_angle) and (prim_dir < self._max_angle)
        on_limits = (prim_dir == self._min_angle) or (prim_dir == self._max_angle)
        in_range = in_limits or on_limits

        if in_limits:
            halfway = (self._max_angle - self._min_angle) / 2
            if prim_dir <= halfway:
                target = self._min_angle
            else:
                target = self._max_angle
        elif not in_range:
            if prim_dir > self._max_angle:
                target = self._min_angle
            else:
                target = self._max_angle
        else:
            if prim_dir == self._max_angle:
                target = self._min_angle
            else:
                target = self._max_angle


        dur_one_scan = linear_motion_dur(self._min_angle, self._max_angle, prim_slew, prim_accel)

        scan_count = -1 if in_limits else 0
        while True:
            if self._is_dur_based and (time_elapsed >= self._duration - 1e-6):
                break
            elif (not self._is_dur_based) and (scan_count >= self._count):
                break
            
            if self._is_dur_based and (time_elapsed + dur_one_scan) > self._duration:
                time_left = (self._duration - time_elapsed)
                sign = -1 if target < prim_dir else 1
                target = xf_linear_motion_from_time(prim_dir, prim_slew * sign, prim_accel * sign, time_left)


            x, t = linear_motion(prim_dir, target, prim_slew, prim_accel, self._dt)
            y = np.full(len(x), fixed_dir)

            prim_scan_vals.append(x[1:] if len(prim_scan_vals) != 0 else x)
            fixed_scan_vals.append(y[1:] if len(fixed_scan_vals) != 0 else y)
            times.append(times[-1][-1] + t[1:] if len(times) != 0 else t)


            prim_dir = target
            target = self._max_angle if target == self._min_angle else self._min_angle

            time_elapsed += t[-1]
            scan_count += 1
        

        prim_scan_vals = np.concatenate(prim_scan_vals)
        fixed_scan_vals = np.concatenate(fixed_scan_vals)
        times = np.concatenate(times)

        return (prim_scan_vals if is_az_scan else fixed_scan_vals, fixed_scan_vals if is_az_scan else prim_scan_vals, times), (scan_count, time_elapsed)