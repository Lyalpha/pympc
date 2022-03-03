import datetime
import gzip
import logging
import os
import shutil
import tempfile
import urllib.request
from itertools import repeat
from math import pi
from multiprocessing import Pool
from time import time

import astropy.units as u
import ephem
import numpy as np
import pandas as pd
from astropy.time import Time

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s  %(levelname)-10s %(processName)s - %(message)s"
)
logger = logging.getLogger(__name__)

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
        fd, temp_filepath = tempfile.mkstemp(suffix=os.path.splitext(filename)[1], dir=cat_dir)
        response = urllib.request.urlopen(url)
        with open(temp_filepath, "wb") as f:
            f.write(gzip.decompress(response.read()))
        filepath = os.path.join(os.path.dirname(temp_filepath), filename)
        shutil.move(temp_filepath, filepath)
        return filepath

    cats_to_process = [("mpcorb", 0)]
    for include, additional_cat in [(include_nea, ("nea", 1)), (include_comets, ("comets", 2))]:
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

    xephem_csv_filepath = _generate_mpcorb_xephem(
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


def _generate_mpcorb_xephem(mpcorb_filepath, nea_filepath=None, comet_filepath=None):
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
        comet_json.rename(columns={"Designation_and_name": "Principal_desig"}, inplace=True)
        # semi-major axis
        comet_json["a"] = comet_json["Perihelion_dist"] / (1 - comet_json["e"])
        # mean daily motion
        comet_json["n"] = 360.0 / (365.25689 * comet_json["a"] ** 1.5)

        # epoch (as JD)
        def to_jd(row, year, month, day):
            return (
                Time("{}-{:02}-{:02}".format(*map(int, (row[year], row[month], row[day]))))
                + (row[day] % 1) * u.day
            ).jd

        comet_json["Epoch"] = comet_json.apply(
            to_jd, axis=1, args=("Epoch_year", "Epoch_month", "Epoch_day")
        )
        # mean anomoly
        comet_json["PeriEpoch"] = comet_json.apply(
            to_jd, axis=1, args=("Year_of_perihelion", "Month_of_perihelion", "Day_of_perihelion")
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
        mpcorb_json = mpcorb_json[~mpcorb_json["Principal_desig"].isin(nea_json["Principal_desig"])]
        # append the catalogue rows onto this cut-down mpcorb
        mpcorb_json = pd.concat([mpcorb_json, nea_json], sort=False)

    logger.info("creating xephem format database from mpcorb catalogue")
    # Where we don't have a "Name" for the object, we use the "Prinicpal_desig"
    mpcorb_json["Name"] = mpcorb_json["Name"].mask(pd.isnull, mpcorb_json["Principal_desig"])
    # write a minimal version of the catalogue in xephem format - column order is important
    # eccentric orbits (eccentricity < 1)
    xephem_db_e = mpcorb_json.loc[
        mpcorb_json.e < 1, ["Name", "i", "Node", "Peri", "a", "n", "e", "M", "Epoch", "H", "G"]
    ].copy()
    xephem_db_e.insert(1, "type", "e")
    xephem_db_e.insert(10, "relative_epoch", 2000)
    xephem_db_e.loc[:, "Epoch"] = Time(xephem_db_e.Epoch, format="jd").decimalyear
    # hyperbolic orbits (eccentricity > 1)
    xephem_db_h = mpcorb_json.loc[
        mpcorb_json.e > 1,
        ["Name", "PeriEpoch", "i", "Node", "Peri", "e", "Perihelion_dist", "H", "G"],
    ].copy()
    xephem_db_h.insert(1, "type", "h")
    xephem_db_h.insert(8, "relative_epoch", 2000)
    xephem_db_h.loc[:, "PeriEpoch"] = Time(xephem_db_h.PeriEpoch, format="jd").decimalyear
    # parabolic orbits (eccentricity = 1)
    xephem_db_p = mpcorb_json.loc[
        mpcorb_json.e == 1, ["Name", "PeriEpoch", "i", "Peri", "Perihelion_dist", "Node", "H", "G"]
    ].copy()
    xephem_db_p.insert(1, "type", "p")
    xephem_db_p.insert(7, "relative_epoch", 2000)
    xephem_db_p.loc[:, "PeriEpoch"] = Time(xephem_db_p.PeriEpoch, format="jd").decimalyear

    logger.info("writing mpcorb xephem database")
    xephem_csv_path = os.path.join(os.path.dirname(mpcorb_filepath), MPCORB_XEPHEM)
    for xephem_db, mode in zip((xephem_db_e, xephem_db_h, xephem_db_p), ("w", "a", "a")):
        xephem_db.to_csv(xephem_csv_path, header=False, index=False, float_format="%.8f", mode=mode)
    logger.info("mpcorb xephem csv database saved to {}".format(xephem_csv_path))


def minor_planet_check(
    ra, dec, epoch, search_radius, xephem_filepath=None, max_mag=None, chunk_size=2e4, quiet=False
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
        defaults to searching for the mpcorb xephem db in `/tmp/`
    max_mag : float, optional
        maximum magnitude of minor planet matches to return.
    chunk_size : int, optional
        the chunk size for multiprocessing of the search. avoid setting too low
        (<<1e4) to avoid large setup time costs. set to 0 to disable multiprocessing.
    quiet : bool, optional
        whether to display informational logging

    Returns
    -------
    results : list of length 4-tuples
        a list of matching minor planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body,
        xephem db-formatted string of matched body)
    """
    if quiet:
        logger.setLevel(logging.WARNING)
    coo = []
    for c, name in zip((ra, dec), ("ra", "dec")):
        if isinstance(c, (int, float)):
            coo.append(c * pi / 180.0)
        else:
            try:
                coo.append(c.to(u.radian).value)
            except (u.UnitConversionError, AttributeError):
                logger.error("could not convert {} {} to radians".format(name, c))
                raise

    if isinstance(epoch, (int, float)):
        epoch = Time(epoch, format="mjd")
    elif isinstance(epoch, datetime.datetime):
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
            logger.error("could not convert search_radius {} to arcseconds".format(search_radius))
            raise

    return _minor_planet_check(
        coo[0], coo[1], decimalyear, search_radius, xephem_filepath, max_mag, c=chunk_size
    )


def _minor_planet_check(ra, dec, epoch, search_radius, xephem_filepath=None, max_mag=None, c=2e4):
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
        defaults to searching for the mpcorb xephem db in `/tmp/`
    max_mag : float, optional
        maximum magnitude of minor planet matches to return.
    c : int, optional
        chunk size when multiprocessing. set to 0 to disable multiprocessing.

    Returns
    -------
    results : list of length 4-tuples
        a list of matching minor planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body,
        xephem db-formatted string of matched body)
    """
    if xephem_filepath is None:
        xephem_filepath = os.path.join(tempfile.gettempdir(), MPCORB_XEPHEM)
    try:
        with open(xephem_filepath, "r") as f:
            xephem_db = f.readlines()
    except FileNotFoundError:
        logger.exception(
            "xephem db csv file not found at {}. run pympc.update_catalogue() "
            "if necessary.".format(xephem_filepath)
        )
        raise
    logger.info(
        "searching for minor planets within {:.2f} arcsec of ra, dec = {:.5f}, {:.5f} at MJD = {:.5f}".format(
            search_radius, ra, dec, Time(epoch, format="decimalyear").mjd
        )
    )
    date = ephem.date(str(epoch))
    t0 = time()
    if c == 0:
        results = _cone_search_xephem_entries(xephem_db, (ra, dec), date, search_radius, max_mag)
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
                ),
            )
        # flatten our list of lists
        results = [r for result in results for r in result]
    logger.info("search took {:.1f} seconds".format(time() - t0))
    logger.info("found {} matches".format(len(results)))
    return results


def _cone_search_xephem_entries(xephem_db, coo, date, search_radius, max_mag=None):
    """
    performs a cone search around a `ra`, `dec` position at `date` to locate any entries
    in the provided `xephem_db` entries that match within `search_radius`.

    Parameters
    ----------
    xephem_db : list
        a list of xephem db-formatted strings used to cross match position against
    coo : tuple
        (ra, dec) coordinates of search position in radians
    date : emphem.date
        date at which to search
    search_radius : float
        search radius in arcseconds
    max_mag : float, optional
        maximum magnitude of minor planet matches to return.

    Returns
    -------
    results : list of length 4-tuples
        a list of matching minor planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body,
        xephem db-formatted string of matched body)
    """
    results = []
    radtodeg = 180.0 / pi
    for xephem_str in xephem_db:
        mp = ephem.readdb(xephem_str.strip())
        mp.compute(date)
        separation = 3600.0 * radtodeg * (float(ephem.separation((mp.a_ra, mp.a_dec), coo)))
        if separation <= search_radius and mp.mag <= (max_mag or np.inf):
            results.append(
                [
                    (float(mp.a_ra) * radtodeg, float(mp.a_dec) * radtodeg),
                    separation,
                    mp.mag,
                    xephem_str.strip(),
                ]
            )
    return results


if __name__ == "__main__":
    pass
