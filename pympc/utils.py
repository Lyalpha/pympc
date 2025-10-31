import logging
import os
import sys
import tempfile
from typing import Tuple, Union

import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)

MPC_OBSCODES_URL = "https://data.minorplanetcenter.net/api/obscodes"


def _cache_dir() -> str:
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser(r"~\AppData\Local")
    elif sys.platform == "darwin":
        base = os.path.expanduser("~/Library/Caches")
    else:
        base = os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))
    path = os.path.join(base, "pympc")
    os.makedirs(path, exist_ok=True)
    return path


def _cache_path() -> str:
    return os.path.join(_cache_dir(), "obs_codes.npy")


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
    df = df.sort_values("obscode").drop_duplicates(subset="obscode", keep="last").reset_index(drop=True)

    return _df_to_structured_array(df)


def ensure_obs_codes_cached(update: bool = False) -> str:
    """
    Ensure the obscodes cache exists. If update=True, force re-download.
    Returns the cache file path.
    """
    cache = _cache_path()
    if os.path.exists(cache) and not update:
        return cache

    logger.info("fetching and caching MPC observatory codes")
    obscode_arr = _fetch_obscodes()

    fd, tmp_path = tempfile.mkstemp(dir=_cache_dir(), suffix=".npy")
    os.close(fd)
    try:
        np.save(tmp_path, obscode_arr, allow_pickle=False)
        os.replace(tmp_path, cache)
    finally:
        try:
            os.remove(tmp_path)
        except (FileNotFoundError, PermissionError, OSError):
            pass
        logger.info(f"Saved obscodes cache to {cache}")
        return cache


def _load_obs_codes() -> np.ndarray:
    """
    Load the obscodes structured array, fetching/caching if needed.
    """
    cache = ensure_obs_codes_cached(update=False)
    return np.load(cache, allow_pickle=False)


def get_observatory_data(
    observatory: Union[str, int, Tuple[float, float, float]]
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
            raise ValueError("observatory tuple must be (longitude, rho_cos_phi, rho_sin_phi)")
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
            raise ValueError(f"ambiguous observatory name '{observatory}', please use the obscode.")
    if row is not None:
        return float(row["longitude"]), float(row["rhocosphi"]), float(row["rhosinphi"])
    raise ValueError(f"observatory '{observatory}' not found by code or name.\n"
                     f"try running pympc.utils.ensure_obs_codes_cached(update=True) to refresh the "
                     f"cache and fetch the latest data from the MPC.")