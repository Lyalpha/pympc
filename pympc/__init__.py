from importlib.metadata import PackageNotFoundError, version

from .pympc import (
    get_catalogue_status,
    generate_xephem_catalogue,
    minor_planet_check,
    planet_hill_sphere_check,
    update_catalogue,
)
from .utils import add_logging, update_obscode_cache

try:
    __version__ = version("pympc")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = [
    "generate_xephem_catalogue",
    "get_catalogue_status",
    "minor_planet_check",
    "update_catalogue",
    "planet_hill_sphere_check",
    "add_logging",
    "update_obscode_cache",
]
