from .scansession import ScanSession
from .coveragemap import CoverageMap, Overlay
from .coveragemap.utils import merge_maps
from .utils import format_seconds, save_maps_pdf

import numpy as np

class ScanProject:
    """
    Used to bundle multiple ScanSession instances together
    """
    def __init__(self, sessions: list[ScanSession], project_name: str = "Project"):
        if not isinstance(sessions, list):
            raise ValueError("sessions must be a list of ScanSession instances")
        if not all(isinstance(sess, ScanSession) for sess in sessions):
            raise ValueError("all elements in sessions must be ScanSession instances")
        
        self.project_name = project_name
        self.sessions = [sess.copy() for sess in sessions]

        # initialized on run_scans()
        self.coverage_map: CoverageMap = None

    def run_scans(self, frame: str = "ra_dec"):
        """
        Runs all scans in *sessions*

        Parameters
        ----------
        frame  :  str
          the frame for the scan to be run in ("galactic" or "ra_dec"). "ra_dec" is equivalent to "icrs"
        """

        #reset from previous run
        self.sessions = [sess.copy() for sess in self.sessions]
        self.coverage_map: CoverageMap = None

        if frame not in ["ra_dec", "galactic"]:
            raise ValueError("frame must be either \"ra_dec\" or \"galactic\"")

        for sess in self.sessions:
            sess.run_scans(frame)
        
        sess_maps = [sess.coverage_map for sess in self.sessions]
        
        for sess in self.sessions:
            sess.coverage_map.title = sess.start_time.strftime("%Y-%m-%d %I:%M:%S %p %Z") + " Session Coverage"

        self.coverage_map = merge_maps(sess_maps)
        self.coverage_map.title = self.project_name + " - Coverage"

    def save_project_pdf(self, save_loc: str, threshold: int | None = None, overlays:  None | Overlay | list[Overlay] = None, show_titles: bool = True, show_unit_bars: bool = True):
        """
        Saves all information about the Project as a pdf

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
        project_string = (f"{self.project_name}:\n"
                          f"Total Coverage: {self.coverage_map.get_coverage() * 100:.2f}%\n")
        
        for i, sess in enumerate(self.sessions):
            start = sess.start_time
            end = sess._scans[-1]._fin_altaz.obstime.to_datetime(timezone=sess.start_time.tzinfo)

            session_string = (
                f"\tSession {i + 1}:\n"
                f"\t\tStart time: {start.strftime("%Y-%m-%d %I:%M:%S %p %Z")}\n"
                f"\t\tEnd time: {end.strftime("%Y-%m-%d %I:%M:%S %p %Z")}\n"
                f"\t\tDuration of session: {format_seconds((end-start).total_seconds())}\n"
                f"\t\tSky covered: {sess.coverage_map.get_coverage() * 100:.2f}%\n\n"
            )

            project_string += session_string

        maps = []
        for cv_map in [self.coverage_map] + [sess.coverage_map for sess in self.sessions]:
            cv_map = cv_map.copy()
            cv_map.set_display_pref(show_titles, show_unit_bars, threshold)
            maps.append(cv_map)

        save_maps_pdf(maps, save_loc, project_string, overlays)


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

        for sess in self.sessions:
            sess.coverage_map.set_display_pref(show_titles, threshold=threshold)
            sess.coverage_map.display(overlays)

        

        

