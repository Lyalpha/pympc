# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0]

**Minimum Python version requirement is now 3.9.**

- Added support for the [`ASTORB`](https://asteroid.lowell.edu/astorb/) catalogue as the default base asteroid source. The `source` argument to `update_catalogue()` allows selection between `astorb` or `mpcorb`.
  - ASTORB choice does not preclude overlaying comets and NEA bodies, which are still fetched from the MPC.
- New command-line checking interface (`pympc check`), replacing `minor-planet-check`, now supports custom observatory codes.
- Added `pympc catalogue [update|status]` command-line interface to generate and check status of xephem catalogues from the command line.
- Added `show_progress` parameter to `update_catalogue()` to control display of download progress bars (useful for notebook environments).
  - Also added `--no-progress` command-line flag to `pympc catalogue update` to disable download progress bars.
- Added `pympc update-obscode-cache` command to refresh the locally cached MPC observatory-code cache.
- Generated xephem catalogues now follow the naming pattern `xephem_{source}.csv`, where `{source}` is `astorb|mpcorb`.
- Project configuration consolidated into `pyproject.toml` (setup.cfg and requirements.txt removed).
- Sidereal time is now calculated once per call, rather than per topocentric correction.
- Speed up reading and parsing of xephem CSV file.

- **API changes:**
  - The default base asteroid catalogue is now Lowell Observatory's ASTORB instead of MPCORB.
  - The command-line script `minor-planet-check` has been removed.
  - The command-line interface is now a subcommand-based `pympc` CLI (`pympc check ...`, `pympc catalogue ...`).
  - The command-line arguments `--include-comets` and `--include-nea` have been removed and are now the default. Turn them off with `--no-comets` and `--no-nea`, as needed.
  - `minor_planet_check()` and `planet_hill_sphere_check()` now always return an `astropy.table.Table`, even with no matches (in which the table is empty). Previously an empty list was returned for no matches~~~~.
  - Added `max_workers` argument to `minor_planet_check()` to control the number of parallel workers with a default of `4`.
    Previously, multiprocessing (i.e. `chunk_size > 0`) would use `n_workers = n_cpu`.
  - The default location (`cat_dir`) for storing xephem CSV files is now the user's OS-specific cache directory (e.g. `~/.cache/pympc` on Linux).
    - The xephem csv file(s) is/are now stored in an `xephem` subdirectory of `cat_dir`


- **Documentation:**
  - Updated README for new ASTORB support and revised examples.

- **Tests:**
  - Added ASTORB generation and MPC overlay coverage.

## [1.5.0]

- Catch `RuntimeError` raised by [PyEphem (#239)](https://github.com/brandon-rhodes/pyephem/issues/239) when a
  body is near-parabolic and reemit as a warning instead.
- Remove lingering prints in favour of logging.


## [1.4.0]

- Python 3.8 is minimum supported version.
- Switched Observatory code fetching to use MPC's API and removed pre-bundled observatory codes file.
  `observatory` argument can now be passed also as an IAU-recognized observatory name as well as the 3-letter code.
- Using `loguru` for logging and added `pympc.add_logging()` convenience function to enable logging when
  using `pympc` as an application.
- Progress bars displayed when downloading catalogues.


## [1.3.0] - 2025-01-20
- Added ability to match to Major Solar System bodies and perform a Planet Hill sphere check of coordinates.


## [1.2.2] - 2024-08-28
- Catch `TypeError` for Python < 3.12 when trying to user `importlib.resources.files`.


## [1.2.1] - 2024-08-28
- Fix import statement to be relative for `.utils`
- [DOCS] Updated README examples.


## [1.2.0] - 2024-08-25

- Added support for Python 3.12.
- Fix bug where MPCORB data where epoch data were incorrectly being calulated in UTC instead of TT.
- Minimum `ephem` version is now `4.1.0` to ensure correct computation of major body positions, to appear
  in a future release.
- [TESTS] Added tests for topocentric corrections.
- [DOCS] Updated notebooks and README to reflect improved matching accuracy.


## [1.1.1] - 2022-07-26

- Patch fix to correct a mix-up in tags and local version numbers.


## [1.1.0] - 2022-07-26

- A buffer is applied to `search_radius` when cross-matching prior to calculating topocentric positions. This
  is to prevent missing positional matches that are initially outside the search radius, but for which the topocentric
  correction brings them within it.
- Add `pyerfa` to requirements to catch dubious year warnings
- [DOCS] Added notebook to quantify the degradation of matching with reference epoch separation


## [1.0.0] - 2022-03-04

- Changelog added!
- Added topocentric corrections.
- Uses the user's system temporary directory by default for downloaded catalogues
  instead of hardcoded `/tmp`.
- Add option to `cleanup` argument to `update_catalogue()` to remove downloaded
  catalogues after the xephem database is created (default is `True`).
- The filepath of the created xephem database is now returned from `update_catalogue()`.
- Added command line interface via `minor_planet_check`, once package is installed.
- Results from `pympc.minor_planet_check()` are now formatted in an astropy table and
  an empty search returns `None` instead of an empy list.
- Added first tests.
