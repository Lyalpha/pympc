from .pympc import (
    minor_planet_check,
    update_catalogue,
    planet_hill_sphere_check,
)
from .utils import (
    add_logging,
)


def update_obscode_cache() -> None:
    """
    Update the obscodes cache by re-downloading from the MPC.
    """
    from .utils import ensure_obs_codes_cached

    ensure_obs_codes_cached(update=True)


__version__ = "1.4.0"
