import importlib.resources

import numpy as np

try:
    # Python 3.9+
    __ref = importlib.resources.files(__name__) / "data/obs_codes.npy"
    with importlib.resources.as_file(__ref) as f:
        OBS_CODES_ARRAY = np.load(str(f))
except AttributeError:
    # Python 3.6-3.8
    import pkg_resources
    OBS_CODES_ARRAY_PATH = pkg_resources.resource_filename(
        __name__, "data/obs_codes.npy"
    )
    OBS_CODES_ARRAY = np.load(OBS_CODES_ARRAY_PATH)


def get_observatory_data(obs_code):
    """Return the longitude, rho_cos_phi, rho_sin_phi of an observatory"""
    obs_data = OBS_CODES_ARRAY[OBS_CODES_ARRAY["Code"] == obs_code]
    if len(obs_data) != 1:
        raise ValueError(f"Observatory code {obs_code} not uniquely found")
    obs_data = obs_data[0]
    return obs_data[1], obs_data[2], obs_data[3]
