import numpy as np
import pkg_resources

OBS_CODES_ARRAY_PATH = pkg_resources.resource_filename(__name__, "data/obs_codes.npy")


def get_observatory_data(obs_code):
    """Return the longitude, rho_cos_phi, rho_sin_phi of an observatory"""
    obs_code_array = np.load(OBS_CODES_ARRAY_PATH)
    obs_data = obs_code_array[obs_code_array["Code"] == obs_code]
    if len(obs_data) != 1:
        raise ValueError(f"Observatory code {obs_code} not uniquely found")
    obs_data = obs_data[0]
    return obs_data[1], obs_data[2], obs_data[3]
