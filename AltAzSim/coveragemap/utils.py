import healpy as hp
import numpy as np
from .coveragemap import CoverageMap
from astropy.coordinates import SkyCoord

def convert_nside(maps: CoverageMap | list[CoverageMap], new_nside) -> CoverageMap | list[CoverageMap]:
    """
    Converts the given map from its current nside to a new nside. Given instances are preserved.

    Parameters
    ---------
    maps  :  CoverageMap or list of CoverageMap
        the map(s) to be converted to new_nside
    
    new_nside  :  int
        the nside to be converted to
    
    Returns
    -------
    CoverageMap or list of CoverageMap
        copies of the original instances converted to new_nside
    """
    if not isinstance(maps, list):
        maps = [maps]

    converted = []

    for cv_map in maps:
        to_scale = cv_map.map.copy()
        to_scale[to_scale == hp.UNSEEN] = 0

        scaled = hp.ud_grade(to_scale, new_nside)
        scaled[scaled < 1] = hp.UNSEEN

        new_map = cv_map.copy()
        new_map.nside = new_nside
        new_map.map = scaled

        converted.append(new_map)

    return converted[0] if len(converted) == 1 else converted


def merge_maps(maps: list[CoverageMap], convert_nsides: bool = False) -> CoverageMap:
    """
    Merges all maps into one map

    Parameters
    ---------
    maps : list of CoverageMap
        Must contain at least two CoverageMap instances

    convert_nsides  :  bool 
        If False and all maps are not of the same nside, will raise ValueError. Set to True to convert all CoverageMap instances to the nside of the first CoverageMap instance in maps
    
    Returns
    -------
    CoverageMap
        one CoverageMap instance with the merged map
    """

    if len(maps) < 2:
        raise ValueError("Must provide at least two maps")
    if not all(isinstance(map, CoverageMap) for map in maps):
        raise ValueError("all elements in maps must be of instance CoverageMap")
    if (not convert_nsides) and not all(map.nside == maps[0].nside for map in maps):
        raise ValueError("convert_nsides is False. All CoverageMaps must have the same nside")
    
    if convert_nsides:
        for i, map in enumerate(maps):
            if map.nside != maps[0].nside:
                maps[i] = convert_nside(map, maps[0].nside)
    
    maps = [map.copy() for map in maps]
    for map in maps:
        map.map[map.map == hp.UNSEEN] = 0

    merged_map = merged_map = sum(m.map for m in maps)

    merged = CoverageMap()
    merged.nside = maps[0].nside
    merged.map = merged_map
    return merged

def uncover_maps(scan_map: CoverageMap, maps_to_uncover: CoverageMap | list[CoverageMap]) -> CoverageMap | list[CoverageMap]:
    """
    Work with pre-existing healpy maps, such as the Haslam 408MHz, to uncover using your scan coverage map

    Parameters
    ----------
    scan_map  :  CoverageMap
      The scan coverage map you want use to uncover maps in maps_to_uncover

    maps_to_uncover  :  CoverageMap, list of CoverageMap
      The pre-existing maps you want to be uncovered with scan_map

    Returns
    ------
    CoverageMap
      A list (or single elemnt) of the uncovered pre-existing map(s)
    """
    if not isinstance(maps_to_uncover, list):
        maps_to_uncover = [maps_to_uncover]

    for i, elem in enumerate(maps_to_uncover):
        if not isinstance(elem, CoverageMap):
            raise ValueError(f"Element at index {i} is of type {type(elem).__name__}, expected CoverageMap")
        if elem.nside != scan_map.nside:
            raise ValueError(f"Map at index {i} is of nside {elem.nside}, expected {scan_map.nside}")
    
    uncovered_maps = []
    mask = scan_map.map >= 1
    for map_to_uncover in maps_to_uncover:
        uncovered_map = np.full_like(scan_map.map, hp.UNSEEN)
        uncovered_map[mask] = map_to_uncover.map[mask]
    
        uncovered = map_to_uncover.copy()
        uncovered.title = map_to_uncover.title + " - Uncovered"
        uncovered.map = uncovered_map
        
        uncovered_maps.append(uncovered)

    return uncovered_maps[0] if len(uncovered_maps) == 1 else uncovered_maps

def _populate_maps(coords: SkyCoord, radius_of_scan: np.floating, maps: CoverageMap | list[CoverageMap]):
    """
    used internally for populating coverage maps
    """
    if not isinstance(maps, list):
        maps = [maps]
    if not all(isinstance(map, CoverageMap) for map in maps):
        raise ValueError("all elements in maps must be instances of CoverageMap")
    if not all(map.nside == maps[0].nside for map in maps):
        raise ValueError("all CoverageMap instances in maps must have the same nside")
    
    prev_pxls = np.array([], dtype=int)
    is_gal = coords.frame.name == "galactic"
    for coord in coords:
        long = coord.l.deg if is_gal else coord.ra.deg
        lat = coord.b.deg if is_gal else coord.dec.deg

        theta = np.radians(90 - lat)                    # 0 rads in galactic/ra-dec -> pi/2 in mollview
        phi = np.radians(long)                          # 0 rads in galactic/ra-dec -> 0 in mollview
        obs_vec = hp.ang2vec(theta, phi)
        
        scanned_pxls = hp.query_disc(nside=maps[0].nside, vec=obs_vec, radius=radius_of_scan)
        new_pxls = np.setdiff1d(scanned_pxls, prev_pxls, assume_unique=True)
        
        for map in maps:
            np.add.at(map.map, new_pxls, 1)

        prev_pxls = scanned_pxls