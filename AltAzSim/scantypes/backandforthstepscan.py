from ..utils.linearmotion import linear_motion, linear_motion_dur, xf_linear_motion_from_time
from .scanprofile import ScanProfile
from ..telescope import Telescope

import astropy.units as u
from astropy.units import Quantity
import numpy as np

class BackAndForthStepScan(ScanProfile):
    def _init_internal(self, scan_mode, max_az, min_az, az_slew, max_alt, min_alt, 
                 tip_slew, rotate_secondary_pos, step, duration = None, count = None):
        super()._init_internal(duration, count)

        if not isinstance(scan_mode, str):
            raise ValueError("scan_mode must be a string")
        if scan_mode not in ["az", "tip"]:
            raise ValueError("scan_mode must be either 'az' or 'tip'")
        if not max_az.unit.is_equivalent(u.deg):
            raise ValueError("max_az must have angular units")
        if not min_az.unit.is_equivalent(u.deg):
            raise ValueError("min_az must have angular units")
        if not az_slew.unit.is_equivalent(u.deg/u.s):
            raise ValueError("az_slew must have angular velocity units")
        if not max_alt.unit.is_equivalent(u.deg):
            raise ValueError("max_alt must have angular units")
        if not min_alt.unit.is_equivalent(u.deg):
            raise ValueError("min_alt must have angular units")
        if not tip_slew.unit.is_equivalent(u.deg/u.s):
            raise ValueError("tip_slew must have angular velocity units")
        if not isinstance(rotate_secondary_pos, bool):
            raise ValueError("rotate_secondary_pos must be a boolean")
        if not step.unit.is_equivalent(u.deg):
            raise ValueError("step must have angular units")
        # if step.value <= 0:
        #     raise ValueError("step must be greater than 0")
        
        self._scan_mode = scan_mode
        self._max_az = max_az.to_value(u.deg)
        self._min_az = min_az.to_value(u.deg)
        self._az_slew = az_slew.to_value(u.deg/u.s)
        self._max_alt = max_alt.to_value(u.deg)
        self._min_alt = min_alt.to_value(u.deg)
        self._tip_slew = tip_slew.to_value(u.deg/u.s)
        self._rotate_secondary_pos = rotate_secondary_pos
        self._step = step.to_value(u.deg)


    @classmethod
    def from_duration(cls, scan_mode: str, max_az: Quantity, min_az: Quantity, az_slew: Quantity, max_alt: Quantity, min_alt: Quantity, 
                     tip_slew: Quantity, rotate_secondary_pos: bool, step: Quantity, duration: Quantity):
        self = cls.__new__(cls)
        self._init_internal(scan_mode, max_az, min_az, az_slew, max_alt, min_alt, tip_slew, rotate_secondary_pos, step, duration=duration)
        return self
    
    @classmethod
    def from_count(cls, scan_mode: str, max_az: Quantity, min_az: Quantity, az_slew: Quantity, max_alt: Quantity, min_alt: Quantity, 
                     tip_slew: Quantity, rotate_secondary_pos: bool, step: Quantity, count: int):
        self = cls.__new__(cls)
        self._init_internal(scan_mode, max_az, min_az, az_slew, max_alt, min_alt, tip_slew, rotate_secondary_pos, step, count=count)
        return self
    
    def copy(self):
        copy = self.__class__.__new__(self.__class__)
        copy._init_internal(
            self._scan_mode,
            self._max_az * u.deg,
            self._min_az * u.deg,
            self._az_slew * u.deg/u.s,
            self._max_alt * u.deg,
            self._min_alt * u.deg,
            self._tip_slew * u.deg/u.s,
            self._rotate_secondary_pos,
            self._step * u.deg,
            self._duration * u.s if self._duration is not None else None,
            self._count
        )
        return copy
    
    def validate_for_tel(self, telescope: Telescope):
        if self._max_az < telescope._min_az or  self._max_az > telescope._max_az:
            raise ValueError("max_az is out of range for telescope")
        if self._min_az < telescope._min_az or  self._min_az > telescope._max_az:
            raise ValueError("min_az is out of range for telescope")
        if self._az_slew > telescope._az_slew:
            raise ValueError("az_slew is not valid for telescope")
        if self._max_alt < telescope._min_alt or  self._max_alt > telescope._max_alt:
            raise ValueError("max_alt is out of range for telescope")
        if self._min_alt < telescope._min_alt or  self._min_alt > telescope._max_alt:
            raise ValueError("min_alt is out of range for telescope")
        if self._tip_slew > telescope._tip_slew:
            raise ValueError("tip_slew is not valid for telescope")
    
    def scan_info(self):
        return (
            f"azimuth range: {self._max_az} degs to {self._min_az} degs\n"
            f"azimuth slew: {self._az_slew} degs/s\n"
            f"altitude range: {self._max_alt} degs to {self._min_alt} degs\n"
            f"tip slew: {self._tip_slew} degs/s\n"
            f"{"altitude" if self._scan_mode == "az" else "azimuth"} step: {self._step} degs\n"
            f"{"min to max" if self._rotate_secondary_pos else "max to min"} step direction"
        )

    @property
    def scan_mode(self) -> str:
        return self._scan_mode
    
    @scan_mode.setter
    def scan_mode(self, val):
        raise ValueError("scan_mode cannot be changed after initializaiton")
    
    
    @property
    def max_az(self) -> u.Quantity:
        return self._max_az * u.deg

    @max_az.setter
    def max_az(self, val: u.Quantity):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("max_az must have angular units (deg)")
        self._max_az = val.to_value(u.deg)


    @property
    def min_az(self) -> u.Quantity:
        return self._min_az * u.deg

    @min_az.setter
    def min_az(self, val: u.Quantity):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("min_az must have angular units (deg)")
        self._min_az = val.to_value(u.deg)


    @property
    def az_slew(self) -> u.Quantity:
        return self._az_slew * (u.deg / u.s)

    @az_slew.setter
    def az_slew(self, val: u.Quantity):
        if not val.unit.is_equivalent(u.deg / u.s):
            raise ValueError("az_slew must have angular velocity units (deg/s)")
        self._az_slew = val.to_value(u.deg / u.s)


    @property
    def max_alt(self) -> u.Quantity:
        return self._max_alt * u.deg

    @max_alt.setter
    def max_alt(self, val: u.Quantity):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("max_alt must have angular units (deg)")
        self._max_alt = val.to_value(u.deg)


    @property
    def min_alt(self) -> u.Quantity:
        return self._min_alt * u.deg

    @min_alt.setter
    def min_alt(self, val: u.Quantity):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("min_alt must have angular units (deg)")
        self._min_alt = val.to_value(u.deg)


    @property
    def tip_slew(self) -> u.Quantity:
        return self._tip_slew * (u.deg / u.s)

    @tip_slew.setter
    def tip_slew(self, val: u.Quantity):
        if not val.unit.is_equivalent(u.deg / u.s):
            raise ValueError("tip_slew must have angular velocity units (deg/s)")
        self._tip_slew = val.to_value(u.deg / u.s)


    @property
    def step(self) -> u.Quantity:
        return self._step * u.deg

    @step.setter
    def step(self, val: u.Quantity):
        if not val.unit.is_equivalent(u.deg):
            raise ValueError("step must have angular units (deg)")
        self._step = val.to_value(u.deg)


    @property
    def rotate_secondary_pos(self) -> bool:
        return self._rotate_secondary_pos

    @rotate_secondary_pos.setter
    def rotate_secondary_pos(self, val: bool):
        if not isinstance(val, bool):
            raise ValueError("rotate_secondary_pos must be a boolean")
        self._rotate_secondary_pos = val


    def run(self):
        """
        Slews in one primary direction, then steps in the secondary
        """
        is_az_scan = self._scan_mode == "az"

        prim_dir = self._init_az if is_az_scan else self._init_alt
        prim_max = self._max_az if is_az_scan else self._max_alt
        prim_min = self._min_az if is_az_scan else self._min_alt
        prim_slew = self._az_slew if is_az_scan else self._tip_slew
        prim_accel = self._tel._az_accel if is_az_scan else self._tel._tip_accel
        sec_dir =  self._init_alt if is_az_scan else self._init_az
        sec_max = self._max_alt if is_az_scan else self._max_az
        sec_min = self._min_alt if is_az_scan else self._min_az
        sec_slew = self._tip_slew if is_az_scan else self._az_slew
        sec_accel = self._tel._tip_accel if is_az_scan else self._tel._az_accel

        dur_one_prim_scan = linear_motion_dur(prim_min, prim_max, prim_slew, prim_accel)
        dur_one_sec_scan = linear_motion_dur(sec_min, sec_min + self._step, sec_slew, sec_accel)

        time_elapsed = np.float64(0.0)

        prim_scan_vals = []
        sec_scan_vals = []
        times = []

        # align secondary axis if needed
        if not ((sec_dir == sec_min and self._rotate_secondary_pos) or (sec_dir == sec_max and not self._rotate_secondary_pos)):
            sec_target = sec_min if self._rotate_secondary_pos else sec_max

            x, t = linear_motion(sec_dir, sec_target, sec_slew, sec_accel, self._dt)
            y = np.full(len(x), prim_dir)

            sec_scan_vals.append(x)
            prim_scan_vals.append(y)
            times.append(t)
            
            # if self._config["incl_align_time"]:
            time_elapsed += t[-1]
            if self._is_dur_based and time_elapsed > self._duration:
                raise ValueError("scan does not allot enough time to align secondary axis")
            
            sec_dir = sec_target.copy()


        in_limits = (prim_dir > prim_min) and (prim_dir < prim_max)
        on_limits = (prim_dir == prim_min) or (prim_dir == prim_max)
        in_range = in_limits or on_limits

        if in_limits:
            halfway = (prim_max - prim_min) / 2
            if prim_dir <= halfway:
                prim_target = prim_min
            else:
                prim_target = prim_max
        elif not in_range:
            if prim_dir > prim_max:
                prim_target = prim_min
            else:
                prim_target = prim_max
        else:
            if prim_dir == prim_max:
                prim_target = prim_min
            else:
                prim_target = prim_max
        
        scan_count = -1 if in_limits else 0
        slew_back = False
        while True:
            if self._is_dur_based and (time_elapsed >= self._duration):
                break
            elif (not self._is_dur_based) and (scan_count >= self._count):
                break
            
            # step in primary axis
            if self._is_dur_based and (time_elapsed + dur_one_prim_scan) > self._duration:
                time_left = self._duration - time_elapsed
                sign = -1 if prim_target < prim_dir else 1
                prim_target = xf_linear_motion_from_time(prim_dir, prim_slew * sign, prim_accel * sign, time_left)

            x, t = linear_motion(prim_dir, prim_target, prim_slew, prim_accel, self._dt)
            y = np.full(len(x), sec_dir)

            prim_scan_vals.append(x[1:] if len(prim_scan_vals) != 0 else x)
            sec_scan_vals.append(y[1:] if len(sec_scan_vals) != 0 else y)
            times.append(times[-1][-1] + t[1:] if len(times) != 0 else t)

            prim_dir = prim_target.copy()
            prim_target = prim_max if prim_target == prim_min else prim_min
            
            time_elapsed += t[-1]
            scan_count += 1

            # if we need to slew the secondary axis back to the start
            if slew_back:
                slew_back_speed = self._tel._tip_slew if is_az_scan else self._tel._az_slew

                slew_back_target = sec_min if self._rotate_secondary_pos else sec_max
                slew_back_dur = linear_motion_dur(sec_min, sec_max, sec_slew, sec_accel)

                if self._is_dur_based and (time_elapsed + slew_back_dur) > self._duration:
                    time_left = self._scan_duration - time_elapsed
                    sign = -1 if sec_target < sec_dir else 1
                    slew_back_target = xf_linear_motion_from_time(sec_dir, slew_back_speed * sign, sec_accel * sign, time_left)
                
                x, t = linear_motion(sec_dir, slew_back_target, slew_back_speed, sec_accel, self._dt)
                y = np.full(len(x), prim_dir)

                sec_scan_vals.append(x[1:] if len(sec_scan_vals) != 0 else x)
                prim_scan_vals.append(y[1:] if len(prim_scan_vals) != 0 else y)
                times.append(times[-1][-1] + t[1:] if len(times) != 0 else t)

                sec_dir = slew_back_target.copy()
                time_elapsed += t[-1]

                slew_back = False

            # normally step in secondary axis
            sign = 1 if self._rotate_secondary_pos else -1
            sec_target = sec_dir + self._step * sign

            if sec_target > sec_max:
                sec_target = sec_max
                slew_back = True
            elif sec_target < sec_min:
                sec_target = sec_min
                slew_back = True

            if self._is_dur_based and (time_elapsed + dur_one_sec_scan) > self._duration:
                time_left = self._duration - time_elapsed
                sec_target = xf_linear_motion_from_time(sec_dir, sec_slew * sign, sec_accel * sign, time_left)
            x, t = linear_motion(sec_dir, sec_target, sec_slew, sec_accel, self._dt)
            y = np.full(len(x), prim_dir)

            sec_scan_vals.append(x[1:] if len(sec_scan_vals) != 0 else x)
            prim_scan_vals.append(y[1:] if len(prim_scan_vals) != 0 else y)
            times.append(times[-1][-1] + t[1:] if len(times) != 0 else t)

            sec_dir = sec_target.copy()
            time_elapsed += t[-1]

        az = np.concatenate(prim_scan_vals) if is_az_scan else np.concatenate(sec_scan_vals)
        alt = np.concatenate(sec_scan_vals) if is_az_scan else np.concatenate(prim_scan_vals)

        return (az, alt, np.concatenate(times)), (scan_count, time_elapsed)