import os
import sys
import tempfile
from typing import Tuple, Union

import numpy as np
import pandas as pd
import requests
from loguru import logger
from platformdirs import user_cache_dir

# No logging by default incase being used as a library
logger.disable("pympc")

MPC_OBSCODES_URL = "https://data.minorplanetcenter.net/api/obscodes"


def get_pympc_cache_dir() -> str:
    return user_cache_dir(appname="pympc", ensure_exists=True)


def _df_to_structured_array(df: pd.DataFrame) -> np.ndarray:
    dtype = np.dtype(
        [
            ("obscode", "U3"),
            ("name", "U128"),
            ("short_name", "U64"),
            ("longitude", "f8"),
            ("rhocosphi", "f8"),
            ("rhosinphi", "f8"),
        ]
    )
    arr = np.empty(len(df), dtype=dtype)
    if arr.dtype.names is not None:
        for field in arr.dtype.names:
            arr[field] = df[field].to_numpy()
    order = np.argsort(arr["obscode"])
    return arr[order]


def _fetch_obscodes() -> np.ndarray:
    r = requests.get(MPC_OBSCODES_URL, json={}, timeout=30)
    r.raise_for_status()
    data = r.json()
    df = pd.DataFrame.from_dict(data, orient="index")
    keep = [
        "obscode",
        "name",
        "short_name",
        "longitude",
        "rhocosphi",
        "rhosinphi",
    ]
    df = df[keep]

    for col in ("longitude", "rhocosphi", "rhosinphi"):
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0).astype(float)

    for col in ("obscode", "name", "short_name"):
        df[col] = df[col].fillna("").astype(str)

    # Deduplicate by obscode
    df = (
        df.sort_values("obscode")
        .drop_duplicates(subset="obscode", keep="last")
        .reset_index(drop=True)
    )

    return _df_to_structured_array(df)


def ensure_obs_codes_cached(update: bool = False) -> str:
    """
    Ensure the MPC observatory-code data is available locally.

    On the first call (or when ``update=True``), the current list of
    observatory codes is downloaded from the Minor Planet Center API and
    saved as a NumPy structured array in the OS-specific user-cache
    directory (e.g. ``~/.cache/pympc/obs_codes.npy`` on Linux).  Subsequent
    calls skip the download unless ``update=True`` is passed.

    Parameters
    ----------
    update : bool, optional
        Force a fresh download even if a cached file already exists.
        Defaults to ``False``.

    Returns
    -------
    str
        Absolute path to the local cache file.

    Raises
    ------
    requests.HTTPError
        If the HTTP request to the MPC API fails.

    See Also
    --------
    update_obscode_cache : Convenience wrapper that discards the return value.
    get_observatory_data : Uses the cache to resolve a code or name to coordinates.
    """
    obscodes_cache = os.path.join(get_pympc_cache_dir(), "obs_codes.npy")
    if os.path.exists(obscodes_cache) and not update:
        return obscodes_cache

    logger.info("Fetching and caching MPC observatory codes")
    obscode_arr = _fetch_obscodes()

    fd, tmp_path = tempfile.mkstemp(suffix=".npy")
    os.close(fd)
    try:
        np.save(tmp_path, obscode_arr, allow_pickle=False)
        os.replace(tmp_path, obscodes_cache)
    except BaseException:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise
    logger.info(f"Saved obscodes cache to {obscodes_cache}")
    return obscodes_cache


def _load_obs_codes() -> np.ndarray:
    """
    Load the obscodes structured array, fetching/caching if needed.
    """
    cache = ensure_obs_codes_cached(update=False)
    return np.load(cache, allow_pickle=False)


def get_observatory_data(
    observatory: Union[str, int, Tuple[float, float, float]],
) -> Tuple[float, float, float]:
    """
    Resolve an observatory specification to (longitude [deg], rho_cos_phi, rho_sin_phi).

    - If a length-three tuple of floats is passed, it is returned unchanged.
    - If a 3‑char obscode is passed, it is matched by code.
    - Otherwise, it is matched by name or short name (case‑insensitive) in that order.
    """

    # Tuple passthrough
    if isinstance(observatory, tuple):
        if len(observatory) != 3:
            raise ValueError(
                "observatory tuple must be (longitude, rho_cos_phi, rho_sin_phi)"
            )
        return float(observatory[0]), float(observatory[1]), float(observatory[2])

    if isinstance(observatory, int):
        observatory = str(observatory)

    if not isinstance(observatory, str):
        raise ValueError(f"unrecognised format for observatory: {observatory!r}")

    obs = observatory.strip()
    obscode_arr: np.ndarray = _load_obs_codes()
    row = None

    if len(obs) == 3:
        obs_mask = obscode_arr["obscode"] == obs
        if np.any(obs_mask):
            row = obscode_arr[obs_mask][0]

    if row is None:
        obs_lower = obs.lower()
        name_mask = np.char.lower(obscode_arr["name"]) == obs_lower
        short_mask = np.char.lower(obscode_arr["short_name"]) == obs_lower

        if (name_sum := np.sum(name_mask)) == 1:
            row = obscode_arr[name_mask][0]
        elif (short_sum := np.sum(short_mask)) == 1:
            row = obscode_arr[short_mask][0]
        elif (name_sum + short_sum) > 1:
            raise ValueError(
                f"ambiguous observatory name '{observatory}', please use the obscode."
            )
    if row is not None:
        return float(row["longitude"]), float(row["rhocosphi"]), float(row["rhosinphi"])
    raise ValueError(
        f"observatory '{observatory}' not found by code or name.\n"
        f"try running pympc.utils.ensure_obs_codes_cached(update=True) to refresh the "
        f"cache and fetch the latest data from the MPC."
    )


def add_logging(level="INFO", sink=sys.stderr):
    """
    Enable logging for the application.
    """

    logger.enable("pympc")
    logger.remove()
    logger.add(sys.stderr if sink is None else sink, level=level)


def update_obscode_cache() -> None:
    """
    Download the latest MPC observatory codes and update the local cache.

    This is a convenience wrapper around :func:`ensure_obs_codes_cached`
    that forces a fresh download regardless of whether a cached copy already
    exists.  Call this periodically to pick up new or renamed observatories.

    The cache is stored in the OS-specific user-cache directory:

    * **Linux / other** — ``$XDG_CACHE_HOME/pympc/`` (default
      ``~/.cache/pympc/``)
    * **macOS** — ``~/Library/Caches/pympc/``
    * **Windows** — ``%LOCALAPPDATA%\\pympc\\``

    From the command line this is equivalent to running::

        pympc update-obscode-cache

    Raises
    ------
    requests.HTTPError
        If the HTTP request to the MPC API fails.

    See Also
    --------
    ensure_obs_codes_cached : Lower-level function that skips the download
        if an up-to-date cache already exists.
    """
    ensure_obs_codes_cached(update=True)
