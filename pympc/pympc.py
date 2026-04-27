import glob
import gzip
import os
import shutil
import tempfile
import urllib.request
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from itertools import repeat
from time import time
from typing import TypedDict

import astropy.units as u
import ephem
import erfa
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.time import Time
from loguru import logger
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

from .utils import get_observatory_data

# No logging by default incase being used as a library
logger.disable("pympc")

__all__ = [
    "update_catalogue",
    "generate_xephem_catalogue",
    "minor_planet_check",
    "planet_hill_sphere_check",
]


class PlanetProperties(TypedDict):
    """Type definition for planet properties in MAJOR_BODIES."""

    moons: dict[str, float]
    hill_radius: float


# Hill radii (in au) from JPL DE405
MAJOR_BODIES: dict[str, PlanetProperties] = {
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
    "astorb": {
        "url": "https://ftp.lowell.edu/pub/elgb/astorb.dat.gz",
        "filename": "astorb.dat",
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
DEFAULT_CATALOGUE_SOURCE = "astorb"
XEPHEM_FILENAME_TEMPLATE = "xephem_{source}_{datestr}.csv"
LEGACY_MPCORB_XEPHEM = "mpcorb_xephem.csv"
# keep the old symbol for backward compatibility in internal call sites during migration
MPCORB_XEPHEM = LEGACY_MPCORB_XEPHEM
ASTORB_COLUMNS = {
    # https://asteroid.lowell.edu/astorb/
    "Number": (0, 6),
    "Name": (7, 25),
    "H": (42, 47),
    "G": (48, 53),
    "Epoch_year": (106, 110),
    "Epoch_month": (110, 112),
    "Epoch_day": (112, 114),
    "M": (115, 125),
    "Peri": (126, 136),
    "Node": (137, 147),
    "i": (147, 157),
    "e": (158, 168),
    "a": (168, 181),
}
# Shared column contract that must be returned by every catalogue reader.
CATALOGUE_CONTRACT_COLS = [
    "Number",
    "Name",
    "Principal_desig",
    "Epoch",
    "M",
    "Peri",
    "Node",
    "i",
    "e",
    "a",
    "n",
    "H",
    "G",
]
DAY_IN_YEAR = 365.25689
RADTODEG = 180.0 / np.pi
DEGTORAD = 1 / RADTODEG
# angle subtended by Earth's equatorial radius at a distance of 1 AU
SOLAR_PARALLAX_ARCSEC = 8.794143
SOLAR_PARALLAX_RAD = (SOLAR_PARALLAX_ARCSEC / 3600.0) * DEGTORAD


def _xephem_filename(source, dt=None):
    """Build the xephem output filename as xephem_{source}_{YYYYMMDD}.csv."""
    if dt is None:
        dt = datetime.now(timezone.utc)
    return XEPHEM_FILENAME_TEMPLATE.format(source=source, datestr=dt.strftime("%Y%m%d"))


def _normalise_catalogue_source(catalogue_source):
    source = str(catalogue_source).strip().lower()
    if source not in ("mpcorb", "astorb"):
        raise ValueError(
            f"unsupported catalogue_source '{catalogue_source}'. Expected one of: 'mpcorb', 'astorb'."
        )
    return source


def _clean_identifier_series(series):
    if series is None:
        return None
    cleaned = series.fillna("").astype(str).str.strip()
    cleaned = cleaned.replace("", np.nan)
    return cleaned.str.lower()


def _normalise_number_series(series):
    if series is None:
        return None
    cleaned = series.fillna("").astype(str).str.extract(r"(\d+)", expand=False)
    return cleaned.replace("", np.nan)


def _build_orbit_merge_keys(df):
    number = _normalise_number_series(df["Number"]) if "Number" in df.columns else None
    principal = (
        _clean_identifier_series(df["Principal_desig"])
        if "Principal_desig" in df.columns
        else None
    )
    name = _clean_identifier_series(df["Name"]) if "Name" in df.columns else None

    keys = pd.Series(np.nan, index=df.index, dtype=object)
    if principal is not None:
        keys = principal.map(lambda x: f"id:{x}" if pd.notna(x) else np.nan)
    if name is not None:
        name_keys = name.map(lambda x: f"id:{x}" if pd.notna(x) else np.nan)
        keys = keys.where(keys.notna(), name_keys)
    if number is not None:
        number_keys = number.map(lambda x: f"num:{x}" if pd.notna(x) else np.nan)
        keys = keys.where(number_keys.isna(), number_keys)
    return keys


def _overlay_orbit_updates(base_df, overlay_df):
    if overlay_df is None or overlay_df.empty:
        return base_df
    base_keys = _build_orbit_merge_keys(base_df)
    overlay_keys = _build_orbit_merge_keys(overlay_df).dropna()
    if overlay_keys.empty:
        return pd.concat([base_df, overlay_df], sort=False, ignore_index=True)
    base_df = base_df.loc[~base_keys.isin(set(overlay_keys))]
    return pd.concat([base_df, overlay_df], sort=False, ignore_index=True)


def _to_numeric_series(series):
    return pd.to_numeric(
        series.astype(str).str.strip().replace("", np.nan), errors="coerce"
    )


def _calendar_to_jd(year, month, day, scale="tt"):
    jd = pd.Series(np.nan, index=year.index, dtype=float)
    valid = year.notna() & month.notna() & day.notna()
    valid &= month.between(1, 12) & day.between(1, 31)
    if valid.any():
        date_str = (
            year.loc[valid].astype(int).astype(str).str.zfill(4)
            + "-"
            + month.loc[valid].astype(int).astype(str).str.zfill(2)
            + "-"
            + day.loc[valid].astype(int).astype(str).str.zfill(2)
        )
        parsed = pd.to_datetime(date_str, format="%Y-%m-%d", errors="coerce")
        parsed_valid = ~pd.isna(parsed)
        if np.any(parsed_valid):
            jd.loc[date_str.index[parsed_valid]] = Time(
                date_str.loc[parsed_valid].astype(str).tolist(),
                format="iso",
                scale=scale,
            ).jd
    invalid_count = int((~valid).sum())
    if invalid_count:
        logger.warning(
            f"Dropped {invalid_count} rows with invalid calendar epoch components."
        )
    return jd


def _coerce_numeric_columns(df, columns):
    for column in columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def _read_astorb_catalogue(astorb_filepath):
    """Parse the Lowell Observatory ASTORB fixed-width ASCII catalogue.

    ASTORB contains numbered and unnumbered asteroids only — no comets.
    Each 266-column record is parsed according to the documented column spans
    (``ASTORB_COLUMNS``), and only the orbital elements required by downstream
    xephem processing are returned.

    The epoch of osculation is documented by Lowell as TDT (Terrestrial Dynamical
    Time; equivalent to TT).

    Parameters
    ----------
    astorb_filepath : str
        Path to the decompressed ``astorb.dat`` file.

    Returns
    -------
    pandas.DataFrame
        One row per asteroid. Columns are the shared catalogue contract:
        ``Number``, ``Name``, ``Principal_desig``,
        ``Epoch`` (Julian Date, TT scale),
        ``M`` (mean anomaly, deg), ``Peri`` (argument of perihelion, deg),
        ``Node`` (longitude of ascending node, deg), ``i`` (inclination, deg),
        ``e`` (eccentricity), ``a`` (semi-major axis, AU),
        ``n`` (mean daily motion, deg/day), ``H``, ``G``.

        Rows with missing values in any orbital element are dropped.
        ``H`` and ``G`` may be NaN for some rows; callers are responsible
        for applying defaults (see ``generate_xephem_catalogue``).

    """
    logger.info("Creating xephem format file from astorb catalogue")
    logger.info("Reading {}".format(astorb_filepath))
    df = pd.read_fwf(
        astorb_filepath,
        colspecs=list(ASTORB_COLUMNS.values()),
        names=list(ASTORB_COLUMNS.keys()),
        dtype=str,
    )
    for column in ("Number", "Name"):
        df[column] = df[column].fillna("").astype(str).str.strip()

    numeric_columns = [
        "H",
        "G",
        "Epoch_year",
        "Epoch_month",
        "Epoch_day",
        "M",
        "Peri",
        "Node",
        "i",
        "e",
        "a",
    ]
    for column in numeric_columns:
        df[column] = _to_numeric_series(df[column])

    df["Epoch"] = _calendar_to_jd(df["Epoch_year"], df["Epoch_month"], df["Epoch_day"])
    df["n"] = 360.0 / (DAY_IN_YEAR * df["a"] ** 1.5)
    df["Principal_desig"] = df["Name"]

    df = df.dropna(subset=["Name", "Epoch", "M", "Peri", "Node", "i", "e", "a", "n"])
    return df.reindex(columns=CATALOGUE_CONTRACT_COLS)


def _read_mpcorb_catalogue(mpcorb_filepath):
    """Read the MPC extended JSON orbit catalogue.

    Returns a DataFrame restricted to the shared catalogue contract —
    the same columns as ``_read_astorb_catalogue``:
    ``Number``, ``Name``, ``Principal_desig``,
    ``Epoch`` (Julian Date, TT scale), ``M``, ``Peri``,
    ``Node``, ``i``, ``e``, ``a``, ``n``, ``H``, ``G``.

    Comet-specific columns (``PeriEpoch``, ``Perihelion_dist``, …) are not
    preserved here; they are added exclusively by the comet-overlay step in
    ``generate_xephem_catalogue`` from the MPC comet JSON.  Any base MPCORB
    entries with e ≥ 1 that lack those columns will be dropped by the guard
    logic in ``generate_xephem_catalogue``.

    The epoch of osculation stored in the JSON is a Julian Date in TT
    (Terrestrial Time), numerically identical to ASTORB's documented TDT.

    Parameters
    ----------
    mpcorb_filepath : str
        Path to the decompressed ``mpcorb_extended.json`` file.

    Returns
    -------
    pandas.DataFrame
        One row per object, columns restricted to ``CATALOGUE_CONTRACT_COLS``.
        ``H`` and ``G`` may be NaN for some rows; callers are responsible
        for applying defaults (see ``generate_xephem_catalogue``).
    """
    logger.info("Reading {}".format(mpcorb_filepath))
    df = pd.read_json(mpcorb_filepath)
    if "Number" not in df.columns:
        df["Number"] = ""
    if "Principal_desig" not in df.columns:
        df["Principal_desig"] = df["Name"]
    return df.reindex(columns=CATALOGUE_CONTRACT_COLS)


def _read_base_catalogue(base_filepath, source):
    if source == "astorb":
        return _read_astorb_catalogue(base_filepath)
    return _read_mpcorb_catalogue(base_filepath)


def _resolve_xephem_filepath(
    xephem_filepath=None, cat_dir=None, catalogue_source=DEFAULT_CATALOGUE_SOURCE
):
    if xephem_filepath is not None:
        return xephem_filepath
    source = _normalise_catalogue_source(catalogue_source)
    xephem_dir = cat_dir or tempfile.gettempdir()
    pattern = os.path.join(
        xephem_dir, XEPHEM_FILENAME_TEMPLATE.format(source=source, datestr="*")
    )
    matches = sorted(glob.glob(pattern))
    if matches:
        return matches[-1]
    raise FileNotFoundError(
        f"xephem CSV file not found for source '{source}' in {xephem_dir}.\n"
        "Run `pympc.update_catalogue()`, or `pympc update-catalogue` from the command line, if necessary."
    )


def _iter_xephem_entries(xephem_filepath):
    """Yield stripped xephem database lines from disk."""
    with open(xephem_filepath, "r") as f:
        for line in f:
            yield line.strip()


def _iter_xephem_chunks(xephem_filepath, chunk_size):
    """Yield fixed-size chunks of stripped xephem database lines from disk."""
    chunk = []
    for line in _iter_xephem_entries(xephem_filepath):
        chunk.append(line)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def update_catalogue(
    include_nea=True,
    include_comets=True,
    cat_dir=None,
    cleanup=True,
    catalogue_source=DEFAULT_CATALOGUE_SOURCE,
    show_progress=True,
):
    """
    Download asteroid/comet orbit elements and save as CSV file readable by xephem.

    Must be run prior to doing any minor planet checking.

    Parameters
    ----------
    include_nea : boolean, optional
        If the selected base catalogue is being downloaded and `include_nea=True`,
        the Minor Planet Center's near earth asteroid catalogue will also be
        downloaded and used to update rows in the base catalogue.
    include_comets : boolean, optional
        See `include_nea`, but for the comet catalogue.
    cat_dir : str, optional
        The directory in which to store downloaded catalogues and the
        formatted xephem databases. If `None` will default to the user's
        tmp directory.
    cleanup : boolean, optional
        If `True` will remove the downloaded catalogues after they are
        processed into the xephem CSV file.
    catalogue_source : str, optional
        Base asteroid catalogue to use. Allowed values are 'astorb' and
        'mpcorb'.
    show_progress : boolean, optional
        Whether to show progress bars during catalogue download(s).
    """

    source = _normalise_catalogue_source(catalogue_source)

    def download_cat(url, filename):
        gz_fd, gz_path = tempfile.mkstemp(suffix=".gz", dir=cat_dir)
        os.close(gz_fd)
        final_fd, temp_filepath = tempfile.mkstemp(
            suffix=os.path.splitext(filename)[1], dir=cat_dir
        )
        os.close(final_fd)

        try:
            logger.info(f"Downloading {filename} from {url}")
            response = urllib.request.urlopen(url)
            total = response.getheader("Content-Length")
            try:
                total = int(total) if total is not None else 0
            except (TypeError, ValueError):
                total = 0

            chunk_size = 16 * 1024

            def _write_stream(out_f, on_chunk=None):
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    out_f.write(chunk)
                    if on_chunk is not None:
                        on_chunk(len(chunk))

            if show_progress:
                _progress_cols = [
                    TextColumn("[bold blue]{task.description}"),
                    BarColumn(bar_width=None),
                    "[progress.percentage]{task.percentage:>3.0f}%",
                    DownloadColumn(),
                    TransferSpeedColumn(),
                    TimeRemainingColumn(),
                ]
                with Progress(*_progress_cols, transient=False) as progress:
                    task = progress.add_task(
                        filename, total=total if total > 0 else None
                    )
                    with open(gz_path, "wb") as gz_f:
                        _write_stream(
                            gz_f, on_chunk=lambda n: progress.advance(task, n)
                        )
            else:
                with open(gz_path, "wb") as gz_f:
                    _write_stream(gz_f)

            logger.info(f"Saved {filename} to temporary gz {gz_path}")

            # decompress gz into the temp filepath
            logger.info(f"Decompressing {filename}")
            with gzip.open(gz_path, "rb") as gz_f, open(temp_filepath, "wb") as out_f:
                shutil.copyfileobj(gz_f, out_f)

            filepath = os.path.join(os.path.dirname(temp_filepath), filename)
            shutil.move(temp_filepath, filepath)
            logger.info(f"Saved to {filepath}")
            return filepath
        except Exception:
            for p in (gz_path, temp_filepath):
                try:
                    if os.path.exists(p):
                        os.remove(p)
                except OSError:
                    pass
            raise

    cats_to_process: list[tuple[str, int]] = [(source, 0)]
    for include, additional_cat in [
        (include_nea, ("nea", 1)),
        (include_comets, ("comets", 2)),
    ]:
        if include:
            cats_to_process.append(additional_cat)

    cat_filepaths = [None, None, None]
    for cat_name, idx in cats_to_process:
        logger.info("Processing {} catalogue".format(cat_name))
        catalogue = CATALOGUES[cat_name]
        cat_filename = catalogue["filename"]
        cat_url = catalogue["url"]
        cat_filepath = download_cat(cat_url, cat_filename)

        cat_filepaths[idx] = cat_filepath

    with warnings.catch_warnings():
        # Catch possible dubious year warnings from erfa
        warnings.filterwarnings("ignore", category=erfa.ErfaWarning)
        xephem_csv_filepath = generate_xephem_catalogue(
            base_filepath=cat_filepaths[0],
            nea_filepath=cat_filepaths[1],
            comet_filepath=cat_filepaths[2],
            catalogue_source=source,
        )
    if cleanup:
        for cat_filepath in cat_filepaths:
            if cat_filepath is not None:
                logger.info(f"Removing {cat_filepath}")
                os.remove(cat_filepath)
    return xephem_csv_filepath


def generate_xephem_catalogue(
    base_filepath,
    nea_filepath=None,
    comet_filepath=None,
    catalogue_source=DEFAULT_CATALOGUE_SOURCE,
):
    source = _normalise_catalogue_source(catalogue_source)
    base_df = _read_base_catalogue(base_filepath, source)

    if comet_filepath:
        logger.info("Reading {}".format(comet_filepath))
        comet_df = pd.read_json(comet_filepath)
        logger.info("Updating orbits for comet objects")
        # for comets, we need to update column names and calculate a few
        # remove non-standard orbits
        comet_df = comet_df[comet_df["Orbit_type"].isin(["C", "P", "D"])]
        # if no epoch specified, set as perihelion
        for column1, column2 in (
            ("Epoch_year", "Year_of_perihelion"),
            ("Epoch_month", "Month_of_perihelion"),
            ("Epoch_day", "Day_of_perihelion"),
        ):
            comet_df.loc[:, column1] = comet_df[column1].fillna(comet_df[column2])
        # set the principal desig column, used to compare to main base catalogue
        comet_df.rename(
            columns={"Designation_and_name": "Principal_desig"}, inplace=True
        )
        # semi-major axis
        comet_df["a"] = comet_df["Perihelion_dist"] / (1 - comet_df["e"])
        # mean daily motion
        comet_df["n"] = 360.0 / (365.25689 * comet_df["a"] ** 1.5)

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

        comet_df["Epoch"] = comet_df.apply(
            to_jd, axis=1, args=("Epoch_year", "Epoch_month", "Epoch_day")
        )
        # mean anomaly
        comet_df["PeriEpoch"] = comet_df.apply(
            to_jd,
            axis=1,
            args=("Year_of_perihelion", "Month_of_perihelion", "Day_of_perihelion"),
        )
        dt_days = comet_df["Epoch"] - comet_df["PeriEpoch"]
        comet_df["M"] = comet_df["n"] * dt_days
        base_df = _overlay_orbit_updates(base_df, comet_df)

    if nea_filepath:
        logger.info("Reading {}".format(nea_filepath))
        nea_df = pd.read_json(nea_filepath)
        logger.info("Updating orbits for nea objects")
        base_df = _overlay_orbit_updates(base_df, nea_df)

    # Enforce numeric columns used in filtering/conversions after overlays.
    base_df = _coerce_numeric_columns(
        base_df,
        [
            "Epoch",
            "M",
            "Peri",
            "Node",
            "i",
            "e",
            "a",
            "n",
            "H",
            "G",
            "PeriEpoch",
            "Perihelion_dist",
        ],
    )

    # Where we don't have a "Name" for the object, use the "Principal_desig"
    base_df["Name"] = base_df["Name"].mask(pd.isnull, base_df.get("Principal_desig"))
    # G=0.15 is the IAU H-G photometric system default for unknown slope
    # parameter.
    base_df["G"] = base_df["G"].fillna(0.15)
    base_df["H"] = base_df["H"].fillna(99.99)

    # write a minimal version of the catalogue in xephem format - column order is important
    # eccentric orbits (eccentricity < 1)
    xephem_e = base_df.loc[base_df.e < 1].reindex(
        columns=["Name", "i", "Node", "Peri", "a", "n", "e", "M", "Epoch", "H", "G"]
    )
    before_drop = len(xephem_e)
    xephem_e = xephem_e.dropna(
        subset=["Name", "i", "Node", "Peri", "a", "n", "e", "M", "Epoch", "H", "G"]
    )
    dropped_e = before_drop - len(xephem_e)
    if dropped_e:
        logger.warning(
            f"Dropped {dropped_e} elliptical entries with incomplete numeric elements."
        )
    xephem_e.insert(1, "type", "e")
    xephem_e.insert(10, "relative_epoch", 2000)
    xephem_e.loc[:, "Epoch"] = Time(
        xephem_e.Epoch, format="jd", scale="tt"
    ).utc.decimalyear

    # hyperbolic orbits (eccentricity > 1)
    # PeriEpoch and Perihelion_dist are only present on comet-type rows; bare
    # asteroid entries with e>1 (very rare) are dropped to avoid malformed output.
    hyp_mask = base_df.e > 1
    if hyp_mask.any():
        hyp_candidates = base_df.loc[hyp_mask]
        hyp_complete = (
            hyp_candidates.loc[
                hyp_candidates["PeriEpoch"].notna()
                & hyp_candidates["Perihelion_dist"].notna()
            ]
            if "PeriEpoch" in hyp_candidates.columns
            else pd.DataFrame(columns=hyp_candidates.columns)
        )
        dropped = hyp_mask.sum() - len(hyp_complete)
        if dropped:
            logger.warning(
                f"Dropped {dropped} hyperbolic entries missing PeriEpoch/Perihelion_dist."
            )
        xephem_h = hyp_complete.reindex(
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
    else:
        xephem_h = pd.DataFrame(
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
    xephem_h.insert(1, "type", "h")
    xephem_h.insert(8, "relative_epoch", 2000)
    if not xephem_h.empty:
        xephem_h.loc[:, "PeriEpoch"] = Time(
            xephem_h.PeriEpoch, format="jd", scale="tt"
        ).utc.decimalyear

    # parabolic orbits (eccentricity = 1)
    # Same guard: only comet-type rows have PeriEpoch/Perihelion_dist.
    par_mask = base_df.e == 1
    if par_mask.any():
        par_candidates = base_df.loc[par_mask]
        par_complete = (
            par_candidates.loc[
                par_candidates["PeriEpoch"].notna()
                & par_candidates["Perihelion_dist"].notna()
            ]
            if "PeriEpoch" in par_candidates.columns
            else pd.DataFrame(columns=par_candidates.columns)
        )
        dropped = par_mask.sum() - len(par_complete)
        if dropped:
            logger.warning(
                f"Dropped {dropped} parabolic entries missing PeriEpoch/Perihelion_dist."
            )
        xephem_p = par_complete.reindex(
            columns=[
                "Name",
                "PeriEpoch",
                "i",
                "Peri",
                "Perihelion_dist",
                "Node",
                "H",
                "G",
            ]
        )
    else:
        xephem_p = pd.DataFrame(
            columns=[
                "Name",
                "PeriEpoch",
                "i",
                "Peri",
                "Perihelion_dist",
                "Node",
                "H",
                "G",
            ]
        )
    xephem_p.insert(1, "type", "p")
    xephem_p.insert(7, "relative_epoch", 2000)
    if not xephem_p.empty:
        xephem_p.loc[:, "PeriEpoch"] = Time(xephem_p.PeriEpoch, format="jd").decimalyear

    logger.info("Writing xephem CSV file")
    xephem_csv_path = os.path.join(
        os.path.dirname(base_filepath), _xephem_filename(source)
    )
    for xephem_db, mode in zip((xephem_e, xephem_h, xephem_p), ("w", "a", "a")):
        xephem_db.to_csv(
            xephem_csv_path, header=False, index=False, float_format="%.8f", mode=mode
        )
    logger.info("xephem CSV file saved to {}".format(xephem_csv_path))
    return xephem_csv_path


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
    chunk_size=1000,
    max_workers=4,
    cat_dir=None,
    catalogue_source=DEFAULT_CATALOGUE_SOURCE,
):
    """
    Perform a minor planet check around a search position

    This is a convenience call to _minor_planet_check(), which actually does
    the work, and allows more flexibility in argument format.

    Parameters
    ----------
    ra : ~astropy.units.Quantity or ~astropy.coordinates.angles.Longitude or float
        RA of search position - if float, assumed to be in degrees
    dec: ~astropy.units.Quantity or ~astropy.coordinates.angles.Latitude or float
        Dec of search position - if float, assumed to be in degrees
    epoch : ~astropy.time.Time or ~datetime.datetime or float
        Epoch at which to search - if float, assumed to be MJD format.
    search_radius : ~astropy.units.Quantity or float
        Radius around which to search the position for matching minor planets. If
        given as a float, assumed to be in arcseconds.
    xephem_filepath : str, optional
        The xephem CSV file to use for calculating minor body positions. If None,
        defaults to searching for the latest dated xephem catalogue for the
        requested `catalogue_source` in `cat_dir` or the system tempdir.
    max_mag : float, optional
        Maximum magnitude of minor planet matches to return.
    include_minor_bodies : boolean, optional
        Whether to include minor bodies in the search. (i.e. Asteroids and comets
        from the MPCORB and NEA catalogues).
    include_major_bodies : boolean, optional
        Whether to include major bodies in the search. (i.e. Planets and their
        major moons).
    observatory : str or int or tuple, optional
        The observatory to use for calculating topocentric corrections to
        the minor body positions. This can be given as the observatory code
        (see notes) or a tuple of (longitude [degrees], rhocosphi, rhosinphi). The
        default returns geometric positions.
    chunk_size : int, optional
        The chunk size for multiprocessing of the search. Avoid setting too low
        (<1000) to avoid large setup time costs. Set to 0 to disable multiprocessing.
    max_workers : int, optional
        Number of worker processes to use when multiprocessing is enabled
        (i.e. ``chunk_size > 0``). Set to 0 to use all available CPUs. Ignored
        when ``chunk_size=0``.
    cat_dir : str, optional
        Directory in which to search for the generated xephem catalogue when
        `xephem_filepath` is not explicitly provided.
    catalogue_source : str, optional
        Base catalogue source whose generated xephem catalogue should be used
        when `xephem_filepath` is not explicitly provided.

    Returns
    -------
    results : astropy.table.Table object
        A table of minor planet matches, with columns:
            * name - the name of the minor planet
            * ra - the RA of the minor planet at the given epoch
            * dec - the Dec of the minor planet at the given epoch
            * mag - the magnitude of the minor planet as given by the orbital catalogue
            * separation - the angular separation between the minor planet and the search position in arcseconds
            * xephem_str - the xephem string for the minor planet (for major bodies, this will be incomplete)


    Notes
    -----
    If passing an observatory code or name, it should match one defined in the Minor Planet Center
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
                logger.error("Could not convert {} {} to radians".format(name, c))
                raise

    if isinstance(epoch, (int, float)):
        epoch = Time(epoch, format="mjd")
    elif isinstance(epoch, datetime):
        epoch = Time(epoch)
    if isinstance(epoch, Time):
        decimalyear = epoch.decimalyear
    else:
        logger.error("Unrecognised format for date {}".format(epoch))
        raise ValueError

    if not isinstance(search_radius, (int, float)):
        try:
            search_radius = search_radius.to(u.arcsec).value
        except (u.UnitConversionError, AttributeError):
            logger.error(
                "Could not convert search_radius {} to arcseconds".format(search_radius)
            )
            raise

    longitude, rhocosphi, rhosinphi = get_observatory_data(observatory)

    results = _minor_planet_check(
        coo[0],
        coo[1],
        decimalyear,
        search_radius,
        xephem_filepath,
        max_mag,
        longitude,
        rhocosphi,
        rhosinphi,
        include_minor_bodies,
        include_major_bodies,
        chunk_size,
        max_workers,
        cat_dir,
        catalogue_source,
    )
    return _to_astropy_table(results)


def _minor_planet_check(
    ra,
    dec,
    epoch,
    search_radius,
    xephem_filepath=None,
    max_mag=None,
    longitude=0.0,
    rhocosphi=0.0,
    rhosinphi=0.0,
    include_minor_bodies=True,
    include_major_bodies=False,
    c=20000,
    max_workers=4,
    cat_dir=None,
    catalogue_source=DEFAULT_CATALOGUE_SOURCE,
):
    """
    Runs a minor planet check with strict format of arguments

    Parameters
    ----------
    ra : float
        RA of search position in radians
    dec : float
        Declination of search position in radians
    epoch : float
        Epoch at which to search in decimal years (e.g. 2019.12345)
    search_radius : float
        Search radius in arcseconds
    xephem_filepath : str, optional
        The xephem CSV file to use for calculating minor body positions. If None,
        defaults to searching for the latest dated xephem catalogue for the
        requested `catalogue_source` in `cat_dir` or the system tempdir.
    max_mag : float, optional
        Maximum magnitude of minor planet matches to return.
    longitude: float
        Longitude of observer in degrees.
    rhocosphi: float
        Parallax constant for cosine of the latitude of the observer.
    rhosinphi: float
        Parallax constant for sine of the latitude of the observer.
    include_minor_bodies : boolean, optional
        Whether to include minor bodies in the search. (i.e. Asteroids and comets
        from the MPCORB and NEA catalogues).
    include_major_bodies : boolean, optional
        Whether to include major bodies in the search. (i.e. Planets and their
        major moons).
    c : int, optional
        Chunk size when multiprocessing. Set to 0 to disable multiprocessing.
    max_workers : int, optional
        Number of worker processes when multiprocessing is enabled (``c > 0``).
        Set to 0 to use all available CPUs (passed as ``None`` to
        ``ProcessPoolExecutor``). Ignored when ``c=0``.
    cat_dir : str, optional
        Directory in which to search for the generated xephem catalogue when
        `xephem_filepath` is not explicitly provided.
    catalogue_source : str, optional
        Base catalogue source whose generated xephem catalogue should be used
        when `xephem_filepath` is not explicitly provided.

    Returns
    -------
    results : list of length 4-tuples
        A list of matching minor planet entries. Each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body,
        xephem db-formatted string of matched body)
    """
    if max_mag is None:
        max_mag = np.inf
    use_topocentric = any((longitude, rhocosphi, rhosinphi))
    if use_topocentric:
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
        f"Searching within {search_radius:.2f} arcsec of ra, dec = {ra:.5f},"
        f" {dec:.5f} rad at MJD = {mjd:.5f}"
    )
    date = ephem.date(str(epoch))
    local_sidereal_time = None
    if use_topocentric:
        local_sidereal_time = (
            Time(date.datetime()).sidereal_time("apparent", longitude=longitude).radian
        )
    t0 = time()
    if include_minor_bodies:
        xephem_filepath = _resolve_xephem_filepath(
            xephem_filepath, cat_dir, catalogue_source
        )
        logger.info(
            f"Searching for minor bodies using xephem catalogue at {xephem_filepath}"
        )

        if c == 0:
            results = _cone_search_xephem_entries(
                _iter_xephem_entries(xephem_filepath),
                (ra, dec),
                date,
                search_radius,
                max_mag,
                longitude,
                rhocosphi,
                rhosinphi,
                search_radius_buffer,
                local_sidereal_time,
            )
        else:
            args_iter = (
                _iter_xephem_chunks(xephem_filepath, c),
                repeat((ra, dec)),
                repeat(date),  # send a picklable string
                repeat(search_radius),
                repeat(max_mag),
                repeat(longitude),
                repeat(rhocosphi),
                repeat(rhosinphi),
                repeat(search_radius_buffer),
                repeat(local_sidereal_time),
            )
            with ProcessPoolExecutor(max_workers=max_workers or None) as executor:
                results = list(executor.map(_cone_search_xephem_entries, *args_iter))
            # flatten our list of lists
            results = [r for result in results for r in result]
    else:
        results = []

    if include_major_bodies:
        logger.info("Searching for major bodies")
        results += _cone_search_major_bodies(
            (ra, dec),
            date,
            search_radius,
            max_mag,
            longitude,
            rhocosphi,
            rhosinphi,
            search_radius_buffer,
            local_sidereal_time,
        )

    logger.info("Search took {:.1f} seconds".format(time() - t0))
    if len(results):
        logger.info("Found {} matches".format(len(results)))
    else:
        logger.info("No matches found")
    return results


def _cone_search_xephem_entries(
    xephem_db,
    coo,
    date,
    search_radius,
    max_mag,
    longitude,
    rhocosphi,
    rhosinphi,
    buffer,
    local_sidereal_time=None,
):
    """
    Performs a cone search for sources.

    Parameters
    ----------
    xephem_db : iterable
        Iterable of xephem db-formatted strings used to cross-match position against.
    coo : tuple
        (ra, dec) coordinates of search position in radians.
    date : emphem.date
        Date at which to search.
    search_radius : float
        Search radius in arcseconds.
    max_mag : float
        Maximum magnitude of minor planet matches to return.
    longitude: float
        Longitude of observer in degrees.
    rhocosphi: float
        Parallax constant for cosine of the latitude of the observer.
    rhosinphi: float
        Parallax constant for sine of the latitude of the observer.
    buffer: float
        Buffer to add to search radius to ensure all matches are found even after
        topocentric corrections are applied.
    local_sidereal_time: float, optional
        Precomputed local apparent sidereal time in radians. If not provided and
        topocentric corrections are requested, this is computed once per call.

    Returns
    -------
    results : list of length 4-tuples
        A list of matching minor planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body,
        xephem db-formatted string of matched body).
    """
    use_topocentric = any((longitude, rhocosphi, rhosinphi))
    if use_topocentric and local_sidereal_time is None:
        local_sidereal_time = (
            Time(date.datetime()).sidereal_time("apparent", longitude=longitude).radian
        )
    search_radius_with_buffer_rad = (search_radius + buffer) * DEGTORAD / 3600.0
    search_radius_rad = search_radius * DEGTORAD / 3600.0
    results = []
    for xephem_str in xephem_db:
        mp = ephem.readdb(xephem_str)
        mp.compute(date)

        # The _cone_search() functionality is directly embedded here to save the overhead
        # of calling another function for every single entry in the xephem catalogue.

        if mp.mag > max_mag:
            continue

        try:
            separation_rad = ephem.separation((mp.a_ra, mp.a_dec), coo)
        except RuntimeError as e:
            warnings.warn(
                f"Could not calculate separation for body {mp.name}: '{str(e)}'",
                RuntimeWarning,
                stacklevel=2,
            )
            continue

        if separation_rad > search_radius_with_buffer_rad:
            continue

        ra, dec = float(mp.a_ra), float(mp.a_dec)
        if use_topocentric:
            ra, dec = equitorial_geocentric_to_topocentric(
                mp.a_ra,
                mp.a_dec,
                mp.earth_distance,
                date.datetime(),
                longitude,
                rhocosphi,
                rhosinphi,
                local_sidereal_time,
            )
            separation_rad = float(ephem.separation((ra, dec), coo))
            if separation_rad > search_radius_rad:
                continue

        separation_arcsec = 3600.0 * RADTODEG * separation_rad

        results.append(
            [
                mp.name,
                ra * RADTODEG,
                dec * RADTODEG,
                separation_arcsec,
                mp.mag,
                xephem_str,
            ]
        )

    return results


def _cone_search_major_bodies(
    coo,
    date,
    search_radius,
    max_mag,
    longitude,
    rhocosphi,
    rhosinphi,
    buffer,
    local_sidereal_time=None,
):
    """
    Performs a cone search for major bodies (i.e. Planets and moons).

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
    rhocosphi: float
        Parallax constant for cosine of the latitude of the observer.
    rhosinphi: float
        Parallax constant for sine of the latitude of the observer.
    buffer: float
        Buffer to add to search radius to ensure all matches are found even after
        topocentric corrections are applied.
    local_sidereal_time: float, optional
        Precomputed local apparent sidereal time in radians used for topocentric
        correction.

    Returns
    -------
    results : list of length 4-tuples
        A list of matching major body entries. each list entry is a tuple of the format
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
                rhocosphi,
                rhosinphi,
                buffer,
                local_sidereal_time,
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
            rhocosphi,
            rhosinphi,
            buffer,
            local_sidereal_time,
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
    Perform a planet Hill sphere check around a search position

    This is a convenience call to _planet_hill_sphere_check(), which
    does the work, with more flexibility in argument format.

    Parameters
    ----------
    ra : ~astropy.units.Quantity or ~astropy.coordinates.angles.Longitude or float
        RA of search position. If given as a float, assumed to be in degrees.
    dec: ~astropy.units.Quantity or ~astropy.coordinates.angles.Latitude or float
        Dec of search position. If given as a float, assumed to be in degrees.
    epoch : ~astropy.time.Time or ~datetime.datetime or float
        epoch at which to search. If given as a float, assumed to be MJD format.
    observatory : str or int or tuple, optional
        The observatory to use for calculating topocentric corrections to
        the minor body positions. This can be given as the observatory code
        (see notes) or a tuple of (longitude [degrees], rhocosphi, rhosinphi). The
        default returns geometric positions.

    Returns
    -------
    results : astropy.table.Table object
        A table of planet matches, with columns:
            * name - the name of the planet
            * ra - the RA of the planet at the given epoch
            * dec - the Dec of the planet at the given epoch
            * mag - the magnitude of the planet as given by the orbital catalogue
            * separation - the angular separation between the planet and the search position in arcseconds

    Notes
    -----
    If passing an observatory code or name, it should match one defined in the Minor Planet Center
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
                logger.error("Could not convert {} {} to radians".format(name, c))
                raise

    if isinstance(epoch, (int, float)):
        epoch = Time(epoch, format="mjd")
    elif isinstance(epoch, datetime):
        epoch = Time(epoch)
    if isinstance(epoch, Time):
        decimalyear = epoch.decimalyear
    else:
        logger.error("Unrecognised format for date {}".format(epoch))
        raise ValueError

    longitude, rhocosphi, rhosinphi = get_observatory_data(observatory)

    results = _planet_hill_sphere_check(
        coo[0],
        coo[1],
        decimalyear,
        longitude,
        rhocosphi,
        rhosinphi,
    )
    return _to_astropy_table(results)


def _planet_hill_sphere_check(
    ra,
    dec,
    epoch,
    longitude=0.0,
    rhocosphi=0.0,
    rhosinphi=0.0,
):
    """
    Performs a planet hill sphere check around a search position.

    Parameters
    ----------
    ra : float
        RA of search position in radians.
    dec : float
        Declination of search position in radians.
    epoch : float
        Epoch at which to search in decimal years (e.g. 2019.12345).
    longitude : float
        Longitude of observer in degrees.
    rhocosphi : float
        Parallax constant for cosine of the latitude of the observer.
    rhosinphi : float
        Parallax constant for sine of the latitude of the observer.

    Returns
    -------
    results : list of length 4-tuples
        A list of matching planet entries. each list entry is a tuple of the format
        ((ra [degrees], dec [degrees]), separation in arcseconds, magnitude of body,
        xephem db-formatted string of matched body).
    """
    date = ephem.date(str(epoch))

    if any([longitude, rhocosphi, rhosinphi]):
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
        hill_sphere_radius = (
            np.arctan(properties["hill_radius"] / planet.earth_distance)
            * RADTODEG
            * 3600
        )

        res = _cone_search(
            None,
            planet,
            (ra, dec),
            date,
            hill_sphere_radius,
            np.inf,
            longitude,
            rhocosphi,
            rhosinphi,
            search_radius_buffer,
        )
        if res is not None:
            results.append(res)
    logger.info("Found {} matches".format(len(results)))
    return results


def _cone_search(
    xephem_str,
    body,
    coo,
    date,
    search_radius,
    max_mag,
    longitude,
    rhocosphi,
    rhosinphi,
    buffer,
    local_sidereal_time=None,
):
    if body.mag > max_mag:
        return None
    use_topocentric = any((longitude, rhocosphi, rhosinphi))
    search_radius_with_buffer_rad = (search_radius + buffer) * DEGTORAD / 3600.0
    search_radius_rad = search_radius * DEGTORAD / 3600.0
    try:
        separation_rad = ephem.separation((body.a_ra, body.a_dec), coo)
    except RuntimeError as e:
        warnings.warn(
            f"Could not calculate separation for body {body.name}: '{str(e)}'",
            RuntimeWarning,
            stacklevel=2,
        )
        return None
    # First match geocentric positions against the buffered search radius
    if separation_rad <= search_radius_with_buffer_rad:
        ra, dec = float(body.a_ra), float(body.a_dec)
        if use_topocentric:
            # Perform a topocentric correction
            ra, dec = equitorial_geocentric_to_topocentric(
                body.a_ra,
                body.a_dec,
                body.earth_distance,
                date.datetime(),
                longitude,
                rhocosphi,
                rhosinphi,
                local_sidereal_time,
            )
            separation_rad = float(ephem.separation((ra, dec), coo))
            # Apply the search radius check again, here without the buffer since we now have
            # topocentric coordinates
            if separation_rad > search_radius_rad:
                return None
        separation_arcsec = 3600.0 * RADTODEG * separation_rad
        return [
            body.name,
            ra * RADTODEG,
            dec * RADTODEG,
            separation_arcsec,
            body.mag,
            xephem_str or body.writedb(),
        ]
    return None


def equitorial_geocentric_to_topocentric(
    ra_geo,
    dec_geo,
    dist_au,
    epoch,
    longitude,
    rhocosphi,
    rhosinphi,
    local_sidereal_time=None,
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
    rhocosphi:
        Parallax constant for cosine of the latitude of the observer.
    rhosinphi:
        Parallax constant for sine of the latitude of the observer.
    local_sidereal_time: float, optional
        Precomputed local apparent sidereal time in radians. When omitted, it is
        computed from ``epoch`` and ``longitude``.

    Returns
    -------
    topo_ra, topo_dec: float
        Topocentric right ascension and declination of the object in radians.

    Notes
    -----
    See https://rdrr.io/github/Susarro/arqastwb/src/R/coordinates.R#sym-geocentric2topocentric
    """

    if local_sidereal_time is None:
        local_sidereal_time = (
            Time(epoch).sidereal_time("apparent", longitude=longitude).radian
        )
    local_hour_angle = local_sidereal_time - ra_geo

    parallax = SOLAR_PARALLAX_RAD / dist_au
    incrra = np.arctan2(
        -(rhocosphi * np.sin(parallax) * np.sin(local_hour_angle)),
        np.cos(dec_geo) - rhocosphi * np.sin(parallax) * np.cos(local_hour_angle),
    )
    ra_topo = ra_geo + incrra
    dec_topo = np.arctan2(
        np.cos(incrra) * (np.sin(dec_geo) - rhosinphi * np.sin(parallax)),
        np.cos(dec_geo) - rhocosphi * np.sin(parallax) * np.cos(local_hour_angle),
    )

    return ra_topo, dec_topo


def _to_astropy_table(results):
    """Convert a list of xephem results to an astropy table"""
    names = ["name", "ra", "dec", "separation", "mag", "xephem_str"]
    # If results is empty list, set to None for compatability with
    # astropy versions <7.1.0
    table = Table(rows=results or None, names=names)
    return table
