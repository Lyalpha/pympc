# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.1] - 2024-08-28
 - Fix bug in package import statement causing error when running `pympc` as an installed package.
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