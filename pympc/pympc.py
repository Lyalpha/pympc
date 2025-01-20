import argparse
import gzip
import importlib.resources
import logging
import os
import shutil
import tempfile
import urllib.request
import warnings
from datetime import datetime
from itertools import repeat
from multiprocessing import Pool
from time import time

import astropy.units as u
import ephem
import erfa
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.time import Time

from .utils import get_observatory_data

logger = logging.getLogger(__name__)

# Hill radii (in au) from JPL DE405
MAJOR_BODIES = {
    "Mercury": {"moons": {}, "hill_radius": 0.0012},
    "Venus": {"moons": {}, "hill_radius": 0.0067},
    "Mars": {"moons": {"Phobos": 11.8, "Deimos": 12.9}, "hill_radius": 0.0066},
    "Jupiter": {
        "moons": {"Io": 5.0, "Europa": 5.3, "Ganymede": 4.6, "Callisto": 5.7},
        "hill_radius": 0.3381,
    },
    "Saturn": {
        "moons": {
            "Mimas": 12.9,
            "Enceladus": 11.7,
            "Tethys": 10.2,
            "Dione": 10.4,
            "Rhea": 9.7,
            "Titan": 8.1,
            "Hyperion": 14.1,
            "Iapetus": 11.1,
        },
        "hill_radius": 0.4120,
    },
    "Uranus": {
        "moons": {
            "Miranda": 16.5,
            "Ariel": 14.3,
            "Umbriel": 15.0,
            "Titania": 13.9,
            "Oberon": 14.1,
        },
        "hill_radius": 0.4464,
    },
    "Neptune": {"moons": {}, "hill_radius": 0.7689},
    "Pluto": {"moons": {}, "hill_radius": 0.0401},
}

CATALOGUES = {
    "mpcorb": {
        "url": "https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz",
        "filename": "mpcorb.json",
    },
    "nea": {
        "url": "https://minorplanetcenter.net/Extended_Files/nea_extended.json.gz",
        "filename": "nea.json",
    },
    "comets": {
        "url": "https://www.minorplanetcenter.net/Extended_Files/cometels.json.gz",
        "filename": "cometels.json",
    },
}
MPCORB_XEPHEM = "mpcorb_xephem.csv"
DAY_IN_YEAR = 365.25689
RADTODEG = 180.0 / np.pi
DEGTORAD = 1 / RADTODEG
# angle subtended by Earth's equatorial radius at a distance of 1 AU
SOLAR_PARALLAX_ARCSEC = 8.794143
SOLAR_PARALLAX_RAD = (SOLAR_PARALLAX_ARCSEC / 3600.0) * DEGTORAD

__ref = importlib.resources.files(__package__) / "data/obs_codes.npy"
with importlib.resources.as_file(__ref) as __ref_file:
    OBS_CODE_ARRAY = np.load(__ref_file)


def update_catalogue(include_nea=True, include_comets=True, cat_dir=None, cleanup=True):
    """
    download MPC asteroid orbit elements and save as csv database readable by xephem

    must be run prior to doing any minor planet checking.

    Parameters
    ----------
    include_nea : boolean, optional
        if the `'mpcorb'` catalogue is being downloaded and `include_nea=True`,
        the Minor Planet Center's near earth asteroid catalogue will also be
        downloaded and used to update rows in mpcorb. (the nea catalogue is
        updated more often and this can help provide better ephemerides.)
    include_comets : boolean, optional
        see `include_nea`, but for the comet catalogue.
    cat_dir : str, optional
        the directory in which to store downloaded catalogues and the
        formatted xephem databases. if `None` will default to the user's
        tmp directory.
    cleanup : boolean, optional
        if `True` will remove the downloaded catalogues after they are
        processed into the xephem database csv file.
    """

    def download_cat(url, filename):
        fd, temp_filepath = tempfile.mkstemp(
            suffix=os.path.splitext(filename)[1], dir=cat_dir
        )
        response = urllib.request.urlopen(url)
        with open(temp_filepath, "wb") as f:
            f.write(gzip.decompress(response.read()))
        filepath = os.path.join(os.path.dirname(temp_filepath), filename)
        shutil.move(temp_filepath, filepath)
        return filepath

    cats_to_process = [("mpcorb", 0)]
    for include, additional_cat in [
        (include_nea, ("nea", 1)),
        (include_comets, ("comets", 2)),
    ]:
        if include:
            cats_to_process.insert(0, additional_cat)

    cat_filepaths = [None, None, None]
    for cat_name, idx in cats_to_process:
        logger.info("processing {} catalogue".format(cat_name))
        catalogue = CATALOGUES[cat_name]
        cat_filename = catalogue["filename"]
        cat_url = catalogue["url"]
        logger.info("downloading from {}".format(cat_url))
        cat_filepath = download_cat(cat_url, cat_filename)
        logger.info("saved as {}".format(cat_filepath))

        cat_filepaths[idx] = cat_filepath

    with warnings.catch_warnings():
        # Catch possible dubious year warnings from erfa
        warnings.filterwarnings("ignore", category=erfa.ErfaWarning)
        xephem_csv_filepath = generate_mpcorb_xephem(
            mpcorb_filepath=cat_filepaths[0],
            nea_filepath=cat_filepaths[1],
            comet_filepath=cat_filepaths[2],
        )
    if cleanup:
        for cat_filepath in cat_filepaths:
            if cat_filepath is not None:
                logging.info(f"removing {cat_filepath}")
                os.remove(cat_filepath)
    return xephem_csv_filepath


def generate_mpcorb_xephem(mpcorb_filepath, nea_filepath=None, comet_filepath=None):
    logger.info("creating xephem format database from mpcorb catalogue")
    logger.info("reading {}".format(mpcorb_filepath))
    mpcorb_json = pd.read_json(mpcorb_filepath)

    if comet_filepath:
        logger.info("reading {}".format(comet_filepath))
        comet_json = pd.read_json(comet_filepath)
        logger.info("updating orbits for comet objects")
        # for comets we need to update column names and calculate a few
        # remove non-standard orbits
        comet_json = comet_json[comet_json["Orbit_type"].isin(["C", "P", "D"])]
        # if no epoch specified, set as perihelion
        for column1, column2 in (
            ("Epoch_year", "Year_of_perihelion"),
            ("Epoch_month", "Month_of_perihelion"),
            ("Epoch_day", "Day_of_perihelion"),
        ):
            comet_json[column1].fillna(comet_json[column2], inplace=True)
        # set the principal deisg column, used to compare to main mpcorb catalogue
        comet_json.rename(
            columns={"Designation_and_name": "Principal_desig"}, inplace=True
        )
        # semi-major axis
        comet_json["a"] = comet_json["Perihelion_dist"] / (1 - comet_json["e"])
        # mean daily motion
        comet_json["n"] = 360.0 / (365.25689 * comet_json["a"] ** 1.5)

        # epoch (as JD)
        def to_jd(row, year, month, day):
            return (
                Time(
                    "{}-{:02}-{:02}".format(
                        *map(int, (row[year], row[month], row[day]))
                    )
                )
                + (row[day] % 1) * u.day
            ).jd

        comet_json["Epoch"] = comet_json.apply(
            to_jd, axis=1, args=("Epoch_year", "Epoch_month", "Epoch_day")
        )
        # mean anomoly
        comet_json["PeriEpoch"] = comet_json.apply(
            to_jd,
            axis=1,
            args=("Year_of_perihelion", "Month_of_perihelion", "Day_of_perihelion"),
        )
        dt_days = comet_json["Epoch"] - comet_json["PeriEpoch"]
        comet_json["M"] = comet_json["n"] * dt_days
        # remove rows in mpcorb for which we have an entry in this catalogue
        mpcorb_json = mpcorb_json[
            ~mpcorb_json["Principal_desig"].isin(comet_json["Principal_desig"])
        ]
        # append the catalogue rows onto this cut-down mpcorb
        mpcorb_json = pd.concat([mpcorb_json, comet_json], sort=False)

    if nea_filepath:
        logger.info("reading {}".format(nea_filepath))
        nea_json = pd.read_json(nea_filepath)
        logger.info("updating orbits for nea objects")
        # remove rows in mpcorb for which we have an entry in this catalogue
        mpcorb_json = mpcorb_json[
            ~mpcorb_json["Principal_desig"].isin(nea_json["Principal_desig"])
        ]
        # append the catalogue rows onto this cut-down mpcorb
        mpcorb_json = pd.concat([mpcorb_json, nea_json], sort=False)

    # Where we don't have a "Name" for the object, we use the "Prinicpal_desig"
    mpcorb_json["Name"] = mpcorb_json["Name"].mask(
        pd.isnull, mpcorb_json["Principal_desig"]
    )
    # write a minimal version of the catalogue in xephem format - column order is important
    # eccentric orbits (eccentricity < 1)
    xephem_db_e = mpcorb_json.loc[mpcorb_json.e < 1].reindex(
        columns=["Name", "i", "Node", "Peri", "a", "n", "e", "M", "Epoch", "H", "G"]
    )
    xephem_db_e.insert(1, "type", "e")
    xephem_db_e.insert(10, "relative_epoch", 2000)
    xephem_db_e.loc[:, "Epoch"] = Time(
        xephem_db_e.Epoch, format="jd", scale="tt"
    ).utc.decimalyear
    # hyperbolic orbits (eccentricity > 1)
    xephem_db_h = mpcorb_json.loc[mpcorb_json.e > 1].reindex(
        columns=[
            "Name",
            "PeriEpoch",
            "i",
            "Node",
            "Peri",
            "e",
            "Perihelion_dist",
            "H",
            "G",
        ]
    )
    xephem_db_h.insert(1, "type", "h")
    xephem_db_h.insert(8, "relative_epoch", 2000)
    xephem_db_h.loc[:, "PeriEpoch"] = Time(
        xephem_db_h.PeriEpoch, format="jd", scale="tt"
    ).utc.decimalyear
    # parabolic orbits (eccentricity = 1)
    xephem_db_p = mpcorb_json.loc[mpcorb_json.e == 1].reindex(
        columns=["Name", "PeriEpoch", "i", "Peri", "Perihelion_dist", "Node", "H", "G"]
    )
    xephem_db_p.insert(1, "type", "p")
    xephem_db_p.insert(7, "relative_epoch", 2000)
    xephem_db_p.loc[:, "PeriEpoch"] = Time(
        xephem_db_p.PeriEpoch, format="jd"
    ).decimalyear

    logger.info("writing mpcorb xephem database")
    xephem_csv_path = os.path.join(os.path.dirname(mpcorb_filepath), MPCORB_XEPHEM)
    for xephem_db, mode in zip(
        (xephem_db_e, xephem_db_h, xephem_db_p), ("w", "a", "a")
    ):
        xephem_db.to_csv(
            xephem_csv_path, header=False, index=False, float_format="%.8f", mode=mode
        )
    logger.info("mpcorb xephem csv database saved to {}".format(xephem_csv_path))
    return xephem_csv_path


def _get_observatory_data(obs_code):
    """Return the longitude, rho_sin_phi, rho_cos_phi of an observatory"""
    obs_data = OBS_CODE_ARRAY[OBS_CODE_ARRAY["Code"] == obs_code]
    if len(obs_data) != 1:
        raise ValueError(f"Observatory code {obs_code} not uniquely found")
    obs_data = obs_data[0]
    return obs_data[1], obs_data[2], obs_data[3]


def minor_planet_check(
    ra,
    dec,
    epoch,
    search_radius,
    xephem_filepath=None,
    max_mag=None,
    include_minor_bodies=True,
    include_major_bodies=True,
    observatory=(0.0, 0.0, 0.0),
    chunk_size=20000,
):
    """
    perform a minor planet check around a search position

    this is a convenience call to _minor_planet_check(), which actually does
    the work, and allows more flexibility in argument format.

    Parameters
    ----------
    ra : ~astropy.units.Quantity or ~astropy.coordinates.angles.Longitude or float
        RA of search position - if float, assumed to be in degrees
    dec: ~astropy.units.Quantity or ~astropy.coordinates.angles.Latitude or float
        Dec of search position - if float, assumed to be in degrees
    epoch : ~astropy.time.Time or ~datetime.datetime or float
        epoch at which to search - if float, assumed to be MJD format.
    search_radius : ~astropy.units.Quantity or float
        radius around which to search the position for matching minor planets - if
        float, assumed to be in arcseconds.
    xephem_filepath : str, optional
        the xephem_db to use for calculating minor body positions. if None,
        defaults to searching for the mpcorb xephem db in system tempdir
    max_mag : float, optional
        maximum magnitude of minor planet matches to return.
    include_minor_bodies : boolean, optional
        whether to include minor bodies in the search. (i.e. asteroids and comets
        from the MPCORB and NEA catalogues).
    include_major_bodies : boolean, optional
        whether to include major bodies in the search. (i.e. planets and their
        major moons).
    observatory : str or int or tuple, optional
        the observatory to use for calculating topocentric corrections to
        the minor body positions. This can be given as the observatory code
        (see notes) or a tuple of (longitude [degrees], rho_cos_phi, rho_sin_phi). The
        default returns geometric positions.
    chunk_size : int, optional
        the chunk size for multiprocessing of the search. avoid setting too low
        (<<1e4) to avoid large setup time costs. set to 0 to disable multiprocessing.

    Returns
    -------
    results : astropy.table.Table object
        a table of minor planet matches, with columns:
            * name - the name of the minor planet
            * ra - the RA of the minor planet at the given epoch
            * dec - the Dec of the minor planet at the given epoch
            * mag - the magnitude of the minor planet as given by the orbital catalogue
            * separation - the angular separation between the minor planet and the search position in arcseconds
            * xephem_str - the xephem string for the minor planet (for major bodies, this will be incomplete)


    Notes
    -----
    If passing an observatory code it should match one defined in the Minor Planet Center
    list of obseratory codes (see https://minorplanetcenter.net/iau/lists/ObsCodes.html).
    """
    coo = []
    for c, name in zip((ra, dec), ("ra", "dec")):
        if isinstance(c, (int, float)):
            coo.append(c * DEGTORAD)
        else:
            try:
                coo.append(c.to(u.radian).value)
            except (u.UnitConversionError, AttributeError):
                logger.error("could not convert {} {} to radians".format(name, c))
                raise

    if isinstance(epoch, (int, float)):
        epoch = Time(epoch, format="mjd")
    elif isinstance(epoch, datetime):
        epoch = Time(epoch)
    if isinstance(epoch, Time):
        decimalyear = epoch.decimalyear
    else:
        logger.error("unrecognised format for date {}".format(epoch))
        raise ValueError

    if not isinstance(search_radius, (int, float)):
        try:
            search_radius = search_radius.to(u.arcsec).value
        except (u.UnitConversionError, AttributeError):
            logger.error(
                "could not convert search_radius {} to arcseconds".format(search_radius)
            )
            raise

    if isinstance(observatory, int):
        observatory = str(observatory)
    if isinstance(observatory, str):
        longitude, rho_cos_phi, rho_sin_phi = get_observatory_data(observatory)
    elif isinstance(observatory, tuple):
        longitude, rho_cos_phi, rho_sin_phi = observatory
    else:
        logger.error(f"unrecognised format for observatory {observatory}")
        raise ValueError

    results = _minor_planet_check(
        coo[0],
        coo[1],
        decimalyear,
        search_radius,
        xephem_filepath,
        max_mag,
        longitude,
        rho_cos_phi,
        rho_sin_phi,
        include_minor_bodies,
        include_major_bodies,
        chunk_size,
    )
    if len(results):
        results = _to_astropy_table(results)
    return results


def _minor_planet_check(
    ra,
    dec,
    epoch,
    search_radius,
    xephem_filepath=None,
    max_mag=None,
    longitude=0.0,
    rho_cos_phi=0.0,
    rho_sin_phi=0.0,
    include_minor_bodies=True,
    include_major_bodies=False,
    c=2e4,
):
    """
    actually runs the minor planet check with strict format of arguments

    Parameters
    ----------
    ra : float
        RA of search position in radians
    dec : float
        Declination of search position in radians
    epoch : float
        epoch at which to search in decimal years (e.g. 2019.12345)
    search_radius : float
        search radius in arcseconds
    xephem_filepath : str, optional
        the xephem_db to use for calculating minor body positions. if None,
        defaults to searching for the mpcorb xephem db in system tempdir
    max_mag : float, optional
        maximum magnitude of minor planet matches to return.
    longitude: float
        Longitude of observer in degrees.
    rho_cos_phi: float
        Parallax constant for cosine of the latitude of the observer.
    rho_sin_phi: float
        Parallax constant for sine of the latitude of the observer.
    include_minor_bodies : boolean, optional
        whether to include minor bodies in the search. (i.e. asteroids and comets
        from the MPCORB and NEA catalogues).
    include_major_bodies : boolean, optional
        whether to include major bodies in the search. (i.e. planets and their
        major moons).
    c : int, optional
        chunk size when multiprocessing. set to 0 to disable multiprocessing.

    Returns
    -------
    results : list of length 4-tuples
        a list of matching minor planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body,
        xephem db-formatted string of matched body)
    """
    if max_mag is None:
        max_mag = np.inf
    if any([longitude, rho_cos_phi, rho_sin_phi]):
        # If we are using topocentric coordinates, there needs to be a small buffer
        # on searching to ensure the geocentric coordinates returned by xephem contain
        # all matches within the search radius, once topocentric corrections are applied.
        # This buffer is roughly the maximal correction for a 1/3 AU distance object. We
        # incur a slight speed penalty for this, but it is necessary.
        search_radius_buffer = SOLAR_PARALLAX_ARCSEC * 3
    else:
        search_radius_buffer = 0

    mjd = Time(epoch, format="decimalyear").mjd
    logger.info(
        f"searching within {search_radius:.2f} arcsec of ra, dec = {ra:.5f},"
        f" {dec:.5f} rad at MJD = {mjd:.5f}"
    )
    date = ephem.date(str(epoch))
    t0 = time()
    if include_minor_bodies:
        logger.info("searching for minor bodies")
        if xephem_filepath is None:
            xephem_filepath = os.path.join(tempfile.gettempdir(), MPCORB_XEPHEM)
        try:
            with open(xephem_filepath, "r") as f:
                xephem_db = [l.strip() for l in f.readlines()]
        except FileNotFoundError:
            logger.exception(
                "xephem db csv file not found at {}. run pympc.update_catalogue() "
                "if necessary.".format(xephem_filepath)
            )
            raise
        if c == 0:
            results = _cone_search_xephem_entries(
                xephem_db,
                (ra, dec),
                date,
                search_radius,
                max_mag,
                longitude,
                rho_cos_phi,
                rho_sin_phi,
                search_radius_buffer,
            )
        else:
            xephem_db_chunks = np.array_split(xephem_db, max(len(xephem_db) // c, 1))
            with Pool() as pool:
                results = pool.starmap(
                    _cone_search_xephem_entries,
                    zip(
                        xephem_db_chunks,
                        repeat((ra, dec)),
                        repeat(date),
                        repeat(search_radius),
                        repeat(max_mag),
                        repeat(longitude),
                        repeat(rho_cos_phi),
                        repeat(rho_sin_phi),
                        repeat(search_radius_buffer),
                    ),
                )
            # flatten our list of lists
            results = [r for result in results for r in result]
    else:
        results = []

    if include_major_bodies:
        logger.info("searching for major bodies")
        results += _cone_search_major_bodies(
            (ra, dec),
            date,
            search_radius,
            max_mag,
            longitude,
            rho_cos_phi,
            rho_sin_phi,
            search_radius_buffer,
        )

    logger.info("search took {:.1f} seconds".format(time() - t0))
    if len(results):
        logger.info("found {} matches".format(len(results)))
    else:
        logger.info("no matches found")
    return results


def _cone_search_xephem_entries(
    xephem_db,
    coo,
    date,
    search_radius,
    max_mag,
    longitude,
    rho_cos_phi,
    rho_sin_phi,
    buffer,
):
    """
    performs a cone search around a position `coo` at epoch `date` to locate any entries
    in the provided `xephem_db` entries that match within `search_radius+buffer`.

    Parameters
    ----------
    xephem_db : list
        a list of xephem db-formatted strings used to cross-match position against
    coo : tuple
        (ra, dec) coordinates of search position in radians
    date : emphem.date
        date at which to search
    search_radius : float
        search radius in arcseconds
    max_mag : float
        maximum magnitude of minor planet matches to return.
    longitude: float
        Longitude of observer in degrees.
    rho_cos_phi: float
        Parallax constant for cosine of the latitude of the observer.
    rho_sin_phi: float
        Parallax constant for sine of the latitude of the observer.
    buffer: float
        Buffer to add to search radius to ensure all matches are found even after
        topocentric corrections are applied.

    Returns
    -------
    results : list of length 4-tuples
        a list of matching minor planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body,
        xephem db-formatted string of matched body)
    """
    results = []
    for xephem_str in xephem_db:
        mp = ephem.readdb(xephem_str)
        mp.compute(date)
        res = _cone_search(
            xephem_str,
            mp,
            coo,
            date,
            search_radius,
            max_mag,
            longitude,
            rho_cos_phi,
            rho_sin_phi,
            buffer,
        )
        if res is not None:
            results.append(res)

    return results


def _cone_search_major_bodies(
    coo,
    date,
    search_radius,
    max_mag,
    longitude,
    rho_cos_phi,
    rho_sin_phi,
    buffer,
):
    """
    performs a cone search around a position `coo` at epoch `date` to locate any major bodies
    that match within `search_radius+buffer`.

    Parameters
    ----------
    coo : tuple
        (ra, dec) coordinates of search position in radians
    date : emphem.date
        date at which to search
    search_radius : float
        search radius in arcseconds
    max_mag : float
        maximum magnitude of matches to return.
    longitude: float
        Longitude of observer in degrees.
    rho_cos_phi: float
        Parallax constant for cosine of the latitude of the observer.
    rho_sin_phi: float
        Parallax constant for sine of the latitude of the observer.
    buffer: float
        Buffer to add to search radius to ensure all matches are found even after
        topocentric corrections are applied.

    Returns
    -------
    results : list of length 4-tuples
        a list of matching major body entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body,
        xephem db-formatted string of matched body)
    """
    results = []
    for planet, properties in MAJOR_BODIES.items():
        planet = getattr(ephem, planet)()
        planet.compute(date)
        for moon, mag in properties["moons"].items():
            moon = getattr(ephem, moon)()
            moon.compute(date)
            moon.earth_distance = planet.earth_distance
            moon.mag = mag
            res = _cone_search(
                None,
                moon,
                coo,
                date,
                search_radius,
                max_mag,
                longitude,
                rho_cos_phi,
                rho_sin_phi,
                buffer,
            )
            if res is not None:
                results.append(res)
        res = _cone_search(
            None,
            planet,
            coo,
            date,
            search_radius,
            max_mag,
            longitude,
            rho_cos_phi,
            rho_sin_phi,
            buffer,
        )
        if res is not None:
            results.append(res)

    return results


def planet_hill_sphere_check(
    ra,
    dec,
    epoch,
    observatory=(0.0, 0.0, 0.0),
):
    """
    perform a planet hill sphere check around a search position

    this is a convenience call to _planet_hill_sphere_check(), which actually does
    the work, and allows more flexibility in argument format.

    Parameters
    ----------
    ra : ~astropy.units.Quantity or ~astropy.coordinates.angles.Longitude or float
        RA of search position - if float, assumed to be in degrees
    dec: ~astropy.units.Quantity or ~astropy.coordinates.angles.Latitude or float
        Dec of search position - if float, assumed to be in degrees
    epoch : ~astropy.time.Time or ~datetime.datetime or float
        epoch at which to search - if float, assumed to be MJD format.
    observatory : str or int or tuple, optional
        the observatory to use for calculating topocentric corrections to
        the minor body positions. This can be given as the observatory code
        (see notes) or a tuple of (longitude [degrees], rho_cos_phi, rho_sin_phi). The
        default returns geometric positions.

    Returns
    -------
    results : astropy.table.Table object
        a table of planet matches, with columns:
            * name - the name of the planet
            * ra - the RA of the planet at the given epoch
            * dec - the Dec of the planet at the given epoch
            * mag - the magnitude of the planet as given by the orbital catalogue
            * separation - the angular separation between the planet and the search position in arcseconds

    Notes
    -----
    If passing an observatory code it should match one defined in the Minor Planet Center
    list of obseratory codes (see https://minorplanetcenter.net/iau/lists/ObsCodes.html).
    """
    coo = []
    for c, name in zip((ra, dec), ("ra", "dec")):
        if isinstance(c, (int, float)):
            coo.append(c * DEGTORAD)
        else:
            try:
                coo.append(c.to(u.radian).value)
            except (u.UnitConversionError, AttributeError):
                logger.error("could not convert {} {} to radians".format(name, c))
                raise

    if isinstance(epoch, (int, float)):
        epoch = Time(epoch, format="mjd")
    elif isinstance(epoch, datetime):
        epoch = Time(epoch)
    if isinstance(epoch, Time):
        decimalyear = epoch.decimalyear
    else:
        logger.error("unrecognised format for date {}".format(epoch))
        raise ValueError

    if isinstance(observatory, int):
        observatory = str(observatory)
    if isinstance(observatory, str):
        longitude, rho_cos_phi, rho_sin_phi = get_observatory_data(observatory)
    elif isinstance(observatory, tuple):
        longitude, rho_cos_phi, rho_sin_phi = observatory
    else:
        logger.error(f"unrecognised format for observatory {observatory}")
        raise ValueError

    results = _planet_hill_sphere_check(
        coo[0],
        coo[1],
        decimalyear,
        longitude,
        rho_cos_phi,
        rho_sin_phi,
    )
    if len(results):
        results = _to_astropy_table(results)
    return results


def _planet_hill_sphere_check(
    ra,
    dec,
    epoch,
    longitude=0.0,
    rho_cos_phi=0.0,
    rho_sin_phi=0.0,
):
    """
    performs a planet hill sphere check around a position `coo` at epoch `date` to locate any major bodies

    Parameters
    ----------
    ra : float
        RA of search position in radians
    dec : float
        Declination of search position in radians
    epoch : float
        epoch at which to search in decimal years (e.g. 2019.12345)
    longitude : float
        Longitude of observer in degrees.
    rho_cos_phi : float
        Parallax constant for cosine of the latitude of the observer.
    rho_sin_phi : float
        Parallax constant for sine of the latitude of the observer.

    Returns
    -------

    """
    date = ephem.date(str(epoch))

    if any([longitude, rho_cos_phi, rho_sin_phi]):
        # Buffer against missing matches due to topocentric corrections
        search_radius_buffer = SOLAR_PARALLAX_ARCSEC * 3
    else:
        search_radius_buffer = 0

    results = []
    for planet, properties in MAJOR_BODIES.items():
        planet = getattr(ephem, planet)()
        planet.compute(date)
        # Calculate hill sphere radius in arcseconds based on value in AU in MAJOR_BODIES dict, and using
        # planet.earth_distance as the distance to the planet in AU
        hill_sphere_radius = np.arctan(
                properties["hill_radius"] / planet.earth_distance
            ) * RADTODEG * 3600

        res = _cone_search(
            None,
            planet,
            (ra, dec),
            date,
            hill_sphere_radius,
            np.inf,
            longitude,
            rho_cos_phi,
            rho_sin_phi,
            search_radius_buffer,
        )
        if res is not None:
            results.append(res)
    logger.info("found {} matches".format(len(results)))
    return results


def _cone_search(
    xephem_str,
    body,
    coo,
    date,
    search_radius,
    max_mag,
    longitude,
    rho_cos_phi,
    rho_sin_phi,
    buffer,
):
    separation = (
        3600.0 * RADTODEG * (float(ephem.separation((body.a_ra, body.a_dec), coo)))
    )
    # First match geocentric positions against the buffered search radius
    if separation <= search_radius + buffer and body.mag <= max_mag:
        ra, dec = float(body.a_ra), float(body.a_dec)
        if any([longitude, rho_cos_phi, rho_sin_phi]):
            # Perform a topocentric correction
            ra, dec = equitorial_geocentric_to_topocentric(
                body.a_ra,
                body.a_dec,
                body.earth_distance,
                date.datetime(),
                longitude,
                rho_cos_phi,
                rho_sin_phi,
            )
            separation = 3600.0 * RADTODEG * (float(ephem.separation((ra, dec), coo)))
            # Apply the search radius check again, here without the buffer since we now have
            # topocentric coordinates
            if separation > search_radius:
                return None
        return [
            body.name,
            ra * RADTODEG,
            dec * RADTODEG,
            separation,
            body.mag,
            xephem_str or body.writedb(),
        ]


def equitorial_geocentric_to_topocentric(
    ra_geo,
    dec_geo,
    dist_au,
    epoch,
    longitude,
    rho_cos_phi,
    rho_sin_phi,
):
    """
    Convert equitorial geocentric coordinates to topocentric coordinates.

    Parameters
    ----------
    ra_geo, dec_geo: float
        Geocentric right ascension and declination of the object in radians.
    dist_au: float
        Distance of the object in AU.
    epoch: datetime.datetime
        Date and time of observation.
    longitude: float
        Longitude of observer in degrees.
    rho_cos_phi:
        Parallax constant for cosine of the latitude of the observer.
    rho_sin_phi:
        Parallax constant for sine of the latitude of the observer.

    Returns
    -------
    topo_ra, topo_dec: float
        Topocentric right ascension and declination of the object in radians.

    Notes
    -----
    See https://rdrr.io/github/Susarro/arqastwb/src/R/coordinates.R#sym-geocentric2topocentric
    """

    local_sidereal_time = (
        Time(epoch).sidereal_time("apparent", longitude=longitude).radian
    )
    local_hour_angle = local_sidereal_time - ra_geo

    parallax = SOLAR_PARALLAX_RAD / dist_au
    incrra = np.arctan2(
        -(rho_cos_phi * np.sin(parallax) * np.sin(local_hour_angle)),
        np.cos(dec_geo) - rho_cos_phi * np.sin(parallax) * np.cos(local_hour_angle),
    )
    ra_topo = ra_geo + incrra
    dec_topo = np.arctan2(
        np.cos(incrra) * (np.sin(dec_geo) - rho_sin_phi * np.sin(parallax)),
        np.cos(dec_geo) - rho_cos_phi * np.sin(parallax) * np.cos(local_hour_angle),
    )

    return ra_topo, dec_topo


def _to_astropy_table(results):
    """Convert a list of xephem results to an astropy table"""
    names = ["name", "ra", "dec", "separation", "mag", "xephem_str"]
    table = Table(rows=results, names=names)
    return table


def _console_script(args=None):
    """Console script to run minor planet checking, optionally downloading catalogues too."""
    parser = argparse.ArgumentParser(
        description="Minor planet checking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "ra",
        type=float,
        help="Right Ascension of position to check in degrees.",
    )

    parser.add_argument(
        "dec",
        type=float,
        help="Declination of position to check in degrees.",
    )
    parser.add_argument(
        "epoch",
        type=str,
        default=None,
        help="Epoch at which to perform the search. Either pass an ISO datetime "
        "string ('YYYY-MM-DDTHH:MM:SS') or Modified Juian Date as a float.",
    )
    parser.add_argument(
        "-r",
        "--radius",
        type=float,
        default=5,
        help="The search radius around ra,dec within which to cross-match minor planets entries.",
    )
    parser.add_argument(
        "--update-mpcorb",
        action="store_true",
        help="Whether to (re-)download the latest MPCORB catalogue. If not given, an existing catalogue "
        "within --cat-dir will be used, and the checking of NEAs or Comets will follow what options "
        "were included when the existing catalogue was downloaded.",
    )
    parser.add_argument(
        "--include-nea",
        action="store_true",
        help="Whether to include Near Earth Asteroids in the catalogue update (this only takes effect when "
        "--update-mpcorb is used).",
    )
    parser.add_argument(
        "--include-comets",
        action="store_true",
        help="Whether to include Comets in the catalogue update (this only takes effect when "
        "--update-mpcorb is used).",
    )
    parser.add_argument(
        "--match-to-major-bodies",
        action="store_true",
        help="Whether to include major bodies (planets and their major moons) in the search.",
    )
    parser.add_argument(
        "--hill-sphere-check",
        action="store_true",
        help="Whether to perform a hill sphere check for major bodies around the search position."
    )
    parser.add_argument(
        "-m",
        "--max-mag",
        type=float,
        default=None,
        help="The maximum magnitude of minor planets to include in the search.",
    )
    parser.add_argument(
        "-o",
        "--observatory",
        type=str,
        default="500",
        help="The three-character observatory code used to define the location of the observer."
        "The default returns geocentric positions.",
    )
    parser.add_argument(
        "-cat-dir",
        type=str,
        default=None,
        help="The directory in which to download catalogues, or fetch an existing catalogue from. "
        "If None, this will default to the current temporary directory. The programme will "
        f"search for a file named {MPCORB_XEPHEM} inside this directory to use.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2e4,
        help="The chunk size for multiprocessing of the search. Avoid setting too low (<<1e4) "
        "to avoid large setup time costs. Set to 0 to disable multiprocessing.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase the verbosity of logging. `-v` for INFO level messages, "
        "`-vv` for DEBUG. Default level is WARNING.",
    )
    args_dict = vars(parser.parse_args(args))

    # Set up logging for command-line usage
    log_levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    log_level = log_levels[min(len(log_levels) - 1, args_dict["verbose"])]
    fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(format=fmt, level=log_level)

    try:
        epoch = float(args_dict["epoch"])
    except ValueError:
        epoch = datetime.fromisoformat(args_dict["epoch"])

    if args_dict["update_mpcorb"]:
        xephem_filepath = update_catalogue(
            args_dict["include_nea"], args_dict["include_comets"], args_dict["cat_dir"]
        )
    else:
        xephem_dir = args_dict["cat_dir"] or tempfile.gettempdir()
        xephem_filepath = os.path.join(xephem_dir, MPCORB_XEPHEM)

    if args_dict["match_to_major_bodies"]:
        print("Major and Minor Planet Check:")
    else:
        print("Minor Planet Check:")
    results = minor_planet_check(
        ra=args_dict["ra"],
        dec=args_dict["dec"],
        epoch=epoch,
        search_radius=args_dict["radius"],
        xephem_filepath=xephem_filepath,
        max_mag=args_dict["max_mag"],
        include_minor_bodies=True,
        include_major_bodies=args_dict["match_to_major_bodies"],
        observatory=args_dict["observatory"],
        chunk_size=args_dict["chunk_size"],
    )
    if len(results):
        del results["xephem_str"]
        print(results.pprint_all())
    else:
        if args_dict["match_to_major_bodies"]:
            print("No major or minor bodies found.")
        else:
            print("No minor planets found.")
    if args_dict["hill_sphere_check"]:
        print("Planet Hill Sphere Check:")
        results_hill = planet_hill_sphere_check(
            ra=args_dict["ra"],
            dec=args_dict["dec"],
            epoch=epoch,
            observatory=args_dict["observatory"],
        )
        if len(results_hill):
            del results_hill["xephem_str"]
            print(results_hill.pprint_all())
        else:
            print("Not found to be inside any planet's hill sphere.")
