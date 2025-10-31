from .pympc import (
    minor_planet_check,
    _minor_planet_check,
    update_catalogue,
    _cone_search_xephem_entries,
    planet_hill_sphere_check,
    _planet_hill_sphere_check,
)


def update_obscode_cache() -> None:
    """
    Update the obscodes cache by re-downloading from the MPC.
    """
    from .utils import ensure_obs_codes_cached

    ensure_obs_codes_cached(update=True)


__version__ = "2.dev0"
