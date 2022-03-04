# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2022-03

 - Changelog added!
 - Uses the user's system temporary directory by default for downloaded catalogues
   instead of hardcoded `/tmp`.
 - Add option to `cleanup` argument to `update_catalogue()` to remove downloaded 
   catalogues after the xephem database is created (default is `True`).
 - The filepath of the created xephem database is now returned from `update_catalogue()`.
 - Added command line interface via `minor_planet_check`, once package is installed.
 - Results from `pympc.minor_planet_check()` are now formatted in an astropy table and 
   an empty search returns `None` instead of an empy list.
 - Added first tests.