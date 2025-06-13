import astropy.units as u
from astropy.units import Quantity
import numpy as np
from abc import ABC, abstractmethod

from ..telescope import Telescope
from ..coveragemap.coveragemap import CoverageMap
from ..config import DT

class ScanProfile(ABC):
    def __init__(self, *args) -> RuntimeError:
        '''
        Use from_duration() or from_count() to initialize
        '''
        
        raise RuntimeError("must use from_duration() or from_count() to initialize")
    
    def _init_internal(self, duration: Quantity | None = None, count: int | None = None):        
        if (duration is None) == (count is None):
            raise ValueError("Specify exactly one of duration or count")
       
        if duration is not None and not duration.unit.is_equivalent(u.s):
            raise ValueError("duration must be given in units time")
        if count is not None and not isinstance(count, int):
            raise ValueError("count must be given as int")
    
        self._duration = duration.to_value(u.s) if duration is not None else None
        self._count = count

        self._is_dur_based = self._duration is not None
        self._dt = DT
       
        # initialized on scan()
        self._tel = None
        self._scan_was_ran = False
        self._coverage_map = None

        self._init_az = None
        self._init_alt = None
        self._init_time = None

        self._fin_az = None
        self._fin_alt = None
        self._fin_time = None
    
    @abstractmethod
    def copy(self):
        pass

    @abstractmethod
    def validate_for_tel(self, telescope: Telescope):
        """
        Ensures valuse used in the scan are valid for the given telescope
        """
        pass

    @abstractmethod
    def scan_info(self) -> str:
        """
        Resturns a description of the scan
        """
        pass
        
    @property
    def duration(self) -> Quantity:
        if not self._is_dur_based and not self._scan_was_ran:
            raise ValueError("No scan has been ran")
        return self._duration * u.s
        
    @duration.setter
    def duration(self, val: Quantity):
        if not self._is_dur_based:
            raise ValueError("Was not initialized as a duration based scan. Cannot mutate")
        if not val.unit.is_equivalent(u.s):
            raise ValueError("duration must have time units (s)")
        self._duration = val.to_value(u.s)
    
    @property
    def count(self) -> int:
        if self._is_dur_based and not self._scan_was_ran:
            raise ValueError("No scan has been ran")
        return self._count

    @count.setter
    def count(self, val: int):
        if self._is_dur_based:
            raise ValueError("Was not initialized as a count based scan. Cannot mutate")
        if not isinstance(val, int):
            raise ValueError("count must be an int")
        self._count = val
    
    @property
    def is_dur_based(self) -> bool:
        return self._is_dur_based
    
    @is_dur_based.setter
    def is_dur_based(self, val):
        raise ValueError("Cannot mutate is_dur_based after initialization")
    
    @property
    def coverage_map(self) -> CoverageMap:
        if not self._scan_was_ran:
            raise ValueError("No scan has been ran")
        return self._coverage_map.copy()
    
    @abstractmethod
    def run(self) -> tuple[tuple[np.floating, np.floating, np.floating], tuple[int, np.floating]]:
        """
        This function is responsible for returning azimuth positions, altitude positions, corresponding times, length of the entire scan, and how many full scans were complete

        return should look like -> ((az, alt, times), (count, dur))
        """
        pass
