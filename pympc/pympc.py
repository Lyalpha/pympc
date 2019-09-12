import datetime
import gzip
import logging
import os
import shutil
import urllib.request
from math import pi
from multiprocessing import Pool
from itertools import repeat
from tempfile import mkstemp
from time import time

import astropy.units as u
import ephem
import numpy as np
import pandas as pd
from astropy.time import Time

logging.basicConfig(level=logging.INFO, format='%(asctime)s  %(levelname)-10s %(processName)s - %(message)s')

MPCORB_EXTENDED_JSON_URL = 'https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz'
DEFAULT_MPCORB_JSON_PATH = '/tmp/mpcorb_extended.json'


def get_mpcorb_json_path():
    return os.environ.get('MPCORB_JSON_PATH', DEFAULT_MPCORB_JSON_PATH)


def get_xephem_csv_path():
    return '{}_xephem.csv'.format(os.path.splitext(get_mpcorb_json_path())[0])


def minor_planet_check(ra, dec, epoch, search_radius, chunk_size=1e4, quiet=False):
    """
    perform a minor planet check around a search position

    this is a convenience call to _minor_planet_check(), which actually does
    the work, and allows more flexibility in argument format.

    Parameters
    ----------
    ra : ~astropy.units.Quantity or ~astropy.coordinates.angles.Longitude or float
        RA of search position, if provided as a float the value is assumed to be
        in degrees
    dec: ~astropy.units.Quantity or ~astropy.coordinates.angles.Latitude or float
        Dec of search position, if provided as a float the value is assumed to be
        in degrees
    epoch : ~astropy.time.Time or ~datetime.datetime or float
        epoch at which to search. If provided as a float the value is assumed to be
        MJD format.
    search_radius : ~astropy.units.Quantity or float
        radius around which to search the position for matching minor planets. If
        provided as a float, the value is assumed to be in arcseconds.
    chunk_size : int
        the chunk size for multiprocessing of the search. avoid setting too low
        (<<1e4) to avoid large setup time costs. set to 0 to disable multiprocessing.
    quiet : bool
        whether to display informational logging

    Returns
    -------
    results : list of length 4-tuples
        a list of matching minor planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body, xephem db-formatted string of matched body)
    """
    if quiet:
        logging.disable(logging.INFO)
    coo = []
    for c, name in zip((ra, dec), ('ra', 'dec')):
        if isinstance(c, (int, float)):
            coo.append(c * pi/180.)
        else:
            try:
                coo.append(c.to(u.radian).value)
            except (u.UnitConversionError, AttributeError):
                logging.error('could not convert {} {} to radians'.format(name, c))
                raise

    if isinstance(epoch, (int, float)):
        epoch = Time(epoch, format='mjd')
    elif isinstance(epoch, datetime.datetime):
        epoch = Time(epoch)
    if isinstance(epoch, Time):
        decimalyear = epoch.decimalyear
    else:
        logging.error('unrecognised format for date {}'.format(epoch))
        raise ValueError

    if not isinstance(search_radius, (int, float)):
        try:
            search_radius = search_radius.to(u.arcsec).value
        except (u.UnitConversionError, AttributeError):
            logging.error('could not convert search_radius {} to arcseconds'.format(search_radius))
            raise

    return _minor_planet_check(coo[0], coo[1], decimalyear, search_radius, c=chunk_size)


def update_catalogue():
    """
    downloads the mpcorb json file and converts it to xephem db format.
    if set, the json file is saved to the path specified in the environment
    variable $MPCORB_JSON_PATH, otherwise it is stored in /tmp/.
    """
    mpcorb_path = get_mpcorb_json_path()
    fd, temp_path = mkstemp(suffix='.json')
    logging.info('downloading mpcorb json catalogue')
    response = urllib.request.urlopen(MPCORB_EXTENDED_JSON_URL)
    with open(temp_path, 'wb') as f:
        f.write(gzip.decompress(response.read()))
    shutil.move(temp_path, mpcorb_path)
    logging.info('mpcorb json catalogue saved as {}'.format(mpcorb_path))

    logging.info('reading mpcorb catalogue')
    mpcorb_json = pd.read_json(mpcorb_path)

    logging.info('creating xephem format database from mpborb catalogue')
    # write a minimal version of the catalogue in xephem format - column order is important
    xephem_db = mpcorb_json[['Name', 'i', 'Node', 'Peri', 'a', 'n', 'e', 'M', 'Epoch', 'H', 'G']].copy()
    xephem_db.insert(1, 'type', 'e')
    xephem_db.insert(10, 'relative_epoch', 2000)
    xephem_db.loc[:, 'Epoch'] = Time(xephem_db.Epoch, format='jd').decimalyear

    logging.info('writing xephem database')
    xephem_csv_path = get_xephem_csv_path()
    xephem_db.to_csv(xephem_csv_path, header=False, index=False)
    logging.info('xephem csv database saved to {}'.format(xephem_csv_path))


def _minor_planet_check(ra, dec, epoch, search_radius, c=1e4):
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
    c : int
        chunk size when multiprocessing. set to 0 to disable multiprocessing.

    Returns
    -------
    results : list of length 4-tuples
        a list of matching minor planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body, xephem db-formatted string of matched body)
    """
    try:
        xephem_db = open(get_xephem_csv_path()).readlines()
    except FileNotFoundError:
        logging.error('xephem csv file not found at {}. set $MPCORB_CAT_PATH as desired and '
                      'run pympc.update_catalogue()'.format(get_xephem_csv_path()))
        return
    date = ephem.date(str(epoch))
    t0 = time()
    if c == 0:
        results = _cone_search_xephem_entries(xephem_db, (ra, dec), date, search_radius)
    else:
        xephem_db_chunks = np.array_split(xephem_db, len(xephem_db)//c)
        with Pool() as pool:
            results = pool.starmap(_cone_search_xephem_entries, zip(xephem_db_chunks, repeat((ra, dec)),
                                                                    repeat(date), repeat(search_radius)))
        # flatten our list of lists
        results = [r for result in results for r in result]
    logging.info('minor planet search took {:.1f} seconds'.format(time() - t0))
    logging.info('found {} matches'.format(len(results)))
    return results


def _cone_search_xephem_entries(xephem_db, coo, date, search_radius=60.):
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

    Returns
    -------
    results : list of length 4-tuples
        a list of matching minor planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body, xephem db-formatted string of matched body)
    """
    results = []
    for xephem_str in xephem_db:
        mp = ephem.readdb(xephem_str.strip())
        mp.compute(date)
        separation = 206264.806*(float(ephem.separation((mp.a_ra, mp.a_dec), coo)))
        if separation < search_radius:
            results.append([(float(mp.a_ra)*180./pi, float(mp.a_dec)*180./pi), separation, mp.mag, xephem_str, ])
    return results


if __name__ == '__main__':
    pass
