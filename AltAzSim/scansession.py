import astropy.units as u
from astropy.units import Quantity
from datetime import datetime
import numpy as np

from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time

from .telescope import Telescope
from .scantypes.scanprofile import ScanProfile
from .coveragemap.overlay import Overlay
from .coveragemap.coveragemap import CoverageMap
from .coveragemap.utils import _populate_maps
from .utils.save_pdf import save_maps_pdf
from .utils.formater import format_seconds

class ScanSession():
    def __init__(self, telescope: Telescope, start_time: datetime, scans: ScanProfile | list[ScanProfile], init_alt: Quantity | None = None, init_az: Quantity | None = None) -> "ScanSession":
        # ensure arguments are of the correct type
        if not isinstance(telescope, Telescope):
            raise ValueError("telescope must be an instance of Telescope")
        if not isinstance(start_time, datetime):
            raise ValueError("start_time must be an instance of datetime")
        if (init_alt is not None) and not init_alt.unit.is_equivalent(u.deg):
            raise ValueError("init_alt must have angular units")
        if (init_az is not None) and not init_az.unit.is_equivalent(u.deg):
            raise ValueError("init_az must have angular units")
        if not isinstance(scans, list):
            scans = [scans]
        if not all(isinstance(scan, ScanProfile) for scan in scans):
            raise TypeError("All elements in scans must be of type ScanProfile")
        
        if init_alt is None: init_alt = telescope._min_alt * u.deg
        if init_az is None: init_az = telescope._min_az * u.deg
        
        # ensure arguments are valid
        if init_alt.to_value(u.deg) > telescope._max_alt or init_alt.to_value(u.deg) < telescope._min_alt:
            raise ValueError("init_alt is not valid for telescope")
        if init_az.to_value(u.deg) > telescope._max_az or init_az.to_value(u.deg) < telescope._min_az:
            raise ValueError("init_az is not valid for telescope")
        
        
        # public members
        self.telescope: Telescope = telescope
        self.start_time: datetime = start_time
        self._scans: list[ScanProfile] = self._init_scans(scans)
        self._init_alt: np.floating = init_alt.to_value(u.deg)
        self._init_az: np.floating = init_az.to_value(u.deg)

        frame = AltAz(obstime=Time(self.start_time), location=self.telescope.location)
        self.init_altaz = SkyCoord(az = init_az, alt= init_alt, frame= frame)

        # initialized on run_scans()
        self.coverage_map: CoverageMap = None
        self.sky_coverage = None

    def copy(self):
        copy = ScanSession(
            self.telescope,
            self.start_time,
            self._init_scans(self._scans),
            self._init_alt * u.deg,
            self._init_az * u.deg
        )

        return copy
    
    @property
    def init_alt(self) -> Quantity:
        return self._init_alt * u.deg
    
    @property
    def init_az(self) -> Quantity:
        return self._init_az * u.deg
    
    @property
    def scans(self) -> ScanProfile | list[ScanProfile]:
        if len(self._scans) > 1:
            return self._scans
        return self._scans[0]
    

    def _init_scans(self, scans: ScanProfile):
        copied_scans = []
        for scan in scans:
            copied = scan.copy()
            copied.validate_for_tel(self.telescope)
            copied_scans.append(copied)
        
        return copied_scans

    def run_scans(self, frame: str = "ra_dec"):
        """
        Runs each scan in *scans*

        Parameters
        ----------
        frame  :  str
          the frame for the scan to be run in ("galactic" or "ra_dec"). "ra_dec" is equivalent to "icrs"

        """
        if frame not in ["ra_dec", "galactic"]:
            raise ValueError("frame must be either \"ra_dec\" or \"galactic\"")
        if frame == "ra_dec":
            frame = 'icrs' #the actual name of the "radec" frame

        # reset from potential previous run_scan call
        self.coverage_map = CoverageMap()
        self._scans = self._init_scans(self._scans)

        cur_az = self._init_az.copy()
        cur_alt = self._init_alt.copy()
        cur_time = Time(self.start_time)

        for i, scan in enumerate(self._scans):
            scan: ScanProfile

            altaz_frame = AltAz(obstime=cur_time, location=self.telescope.location)
            scan._init_alt = cur_alt
            scan._init_az = cur_az
            scan._init_time = cur_time
            scan._scan_was_ran = True
            scan._coverage_map = CoverageMap()
            scan._tel = self.telescope

            (az, alt, times), (count, dur) = scan.run()
            
            # to prevent overlap from where the last scan left off
            if i != 0:
                if times[0] == 0:
                    az = az[1:]
                    alt = alt[1:]
                    times = times[1:]

            if scan._is_dur_based:
                scan._count = count
            else:
                scan._duration = dur
            
            # set the title for the scan
            scan_info = format_seconds(scan._duration) if scan._is_dur_based else f"{scan._count} scans"
            scan._coverage_map.title = f"{cur_time.to_datetime(timezone=self.start_time.tzinfo).strftime("%Y-%m-%d %I:%M:%S %p %Z")} for {scan_info}"

            # transform coords
            obstimes = cur_time + times * u.s
            alt_az_frame = AltAz(obstime=obstimes, location=self.telescope.location)
            coords = SkyCoord(alt=alt * u.deg, az=az * u.deg, frame=alt_az_frame).transform_to(frame)

            _populate_maps(coords, self.telescope._hp_query_radius, [self.coverage_map, scan._coverage_map])

            # set the amount of coverage for the scan
            scan._fract_covered = scan._coverage_map.get_coverage()

            cur_az = az[-1].copy()
            cur_alt = alt[-1].copy()
            cur_time = obstimes[-1].copy()

            # set final conditions for the scan
            altaz_frame = AltAz(obstime=cur_time, location=self.telescope.location)
            scan._fin_altaz = SkyCoord(alt= alt[-1] * u.deg, az= az[-1] * u.deg, frame=altaz_frame)


        # set title for full coverage map
        if len(self._scans) > 1:
            title = "Full Coverage Over all Sessions"
        else:
            title = self._scans[0].coverage_map.title
        
        self.coverage_map.title = title
        self.sky_coverage = self.coverage_map.get_coverage()
        self.fin_altaz = self._scans[-1]._fin_altaz.copy()


    def display_maps(self, threshold: int | None = None, overlays: None | Overlay | list[Overlay] = None, show_titles: bool = True, show_unit_bars: bool = True):
        """
        Displays all maps from the scan session. If using a Jupyter Notebook, maps are displayed inline, otherwise, a popup window opens.

        Parameters
        ---------
        threshold  :  int, optional
          upper bound for the number of passes to be shown. Values greater than threshold apear faded
        
        overlays  :  Overlay or list of Overlay, optional
          point(s) to be overlayed on every map
        
        show_titles  :  bool
          if True, titles are displayed for each coverage map. If False, no title is displayed for every map.

        show_unit_bars  :  bool
          if True, unit bars are displayed for each coverage map. If False, no unit bar is displayed for every map.
        """
        self.coverage_map.set_display_pref(show_titles, show_unit_bars, threshold)
        self.coverage_map.display(overlays)
        if len(self._scans) > 1:
            for scan in self._scans:
                scan._coverage_map.set_display_pref(show_titles, threshold=threshold)
                scan._coverage_map.display(overlays)


    def save_session_pdf(self, save_loc: str, threshold: int | None = None, overlays:  None | Overlay | list[Overlay] = None, show_titles: bool = True, show_unit_bar: bool = True):
        """
        Saves all information about the ScanSession as a pdf

        Parameters
        ---------
        save_loc  :  str
          the save location of pdf
        
        threshold  :  int, optional
          upper bound for the number of passes to be shown. Values greater than threshold apear faded
        
        overlays  :  Overlay or list of Overlay, optional
          point(s) to be overlayed on every map
        
        show_titles  :  bool
          if True, titles are displayed for each coverage map. If False, no title is displayed for every map.

        show_unit_bars  :  bool
          if True, unit bars are displayed for each coverage map. If False, no unit bar is displayed for every map.
        """
        session_string = f"Scan session information:\n\n"
        
        for i, scan in enumerate(self._scans):
            scan_string = (
                f"\tScan {i + 1}:\n"
                f"\t\tStart time: {scan._init_time.to_datetime(timezone=self.start_time.tzinfo).strftime("%Y-%m-%d %I:%M:%S %p %Z")}\n"
                f"\t\tEnd time: {scan._fin_altaz.obstime.to_datetime(timezone=self.start_time.tzinfo).strftime("%Y-%m-%d %I:%M:%S %p %Z")}\n"
                f"\t\tNumber of scans: {scan._count}\n"
                f"\t\tDuration of scan: {format_seconds(scan._duration)}\n"
                f"\t\tSky covered: {scan._fract_covered * 100:.2f}%\n\n"
                f"\t\t{scan.scan_info().replace('\n', '\n\t\t')}\n"
            )

            session_string += scan_string

        if len(self._scans) != 1:
            maps = []
            for cv_map in [self.coverage_map] + [scan._coverage_map for scan in self._scans]:
                cv_map = cv_map.copy()
                cv_map.set_display_pref(show_titles, show_unit_bar, threshold)
                maps.append(cv_map)

        else:
            maps = self.coverage_map.copy()
            maps.set_display_pref(show_titles, show_unit_bar, threshold)
            maps = [maps]

        save_maps_pdf(maps, save_loc, session_string, overlays)

