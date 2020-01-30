pympc
=====

perform checks for the presence of minor bodies at astronomical locations for a given epoch.

#### installation

`pip install .`

or

`python setup.py install`

#### examples

define where we want the orbit catalogue to be downloaded to (see next section)
```
$ export MPCORB_JSON_PATH = /my/file/path/mpcorb_extended.json
```

import the package and (optionally) grab the catalogue - this must be done at least once prior to doing any searches
and can be run to overwrite the catalogue with the latest version from the Minor Planet Center
```
>>> import pympc
>>> pympc.update_catalogue()
```

Note by default it will also download the Near Earth Asteroid catalogue and use orbits from that catalogue where
overlap exists with the main `MPCORB` catalogue, since they are more regularly updated and subject to larger
perturbations. This can be switched off with `pympc.update_catalogue(include_nea=False)`

##### example 1
define our search location, epoch and radius and run the check
```
>>> import astropy.units as u
>>> from astropy.time import Time
>>> ra = 230.028 * u.deg
>>> dec = -11.774 * u.deg
>>> epoch = Time("2019-01-01T00:00")
>>> search_radius = 0.5 * u.arcmin
>>> max_mag = 20
>>> pympc.minor_planet_check(ra, dec, epoch, search_radius, max_mag)
```


##### example 2
here we use float arguments, and the program assumes the units (see `pympc.minor_planet_check()` docstring)
```
>>> ra = 230.028  # assumed degrees
>>> dec = -11.774  # assumed degrees
>>> epoch = 58484.  # assumed MJD
>>> search_radius = 30  # assumed arcseconds
>>> max_mag = 20
>>> pympc.minor_planet_check(ra, dec, epoch, search_radius, max_mag)
```

#### setting the orbit catalogue location
the location to download the MPCORB catalogue is set via the environment variable 
`MPCORB_JSON_PATH` (e.g. `export MPCORB_JSON_PATH = /my/file/path/mpcorb.json`). If this is not found the catalogue
will download to the default location of `/tmp/mpcorb_extended.json`.


#### speed and multiprocessing
the check should take of order a second or two, depending on multiprocessing capabilities.

the private function which actually performs the calculation is `_minor_planet_check()` (note leading underscore).
this can be called to avoid the overhead associated with converting input arguments to `minor_planet_check()`, if
you can provide them directly as required (see `_minor_planet_check()` docstring).

by default the program will calculate positions of bodies in the catalogue multiprocessed. to switch this off set
`chunk_size = 0`, i.e.:

```
>>> pympc.minor_planet_check(ra, dec, epoch, search_radius, chunk_size=0)
```

#### limitations
the orbits are propagated by [xephem](http://www.clearskyinstitute.com/xephem), and that package does not account for
perturbations of the orbits. thus the accuracy of the position is dependant on the time difference between the
epoch of the orbit elements and the epoch at which the search is being performed. 
epoch differences of around a month or less should be fine for typical positional accuracies of ground-based 
observations (~arcsecond), but time differences of several months can accrue to inaccuracies of around ~arcminute for
certain bodies. 

currently the epoch of the orbit elements is visible in the xephem db strings returned by `minor_planet_check()` as a
decimal year format (e.g. ..,2019.317808,..). some diagnostic information and warning when using large time differences
is to be implemented.

the filtering of matches based on magnitude via `max_mag` is limited by the accuracy of the magnitude information in the
database so some buffer should be applied to the desired magnitude cutoff to allow for this.

#### acknowledgments
this package makes use of orbit elements provided by the International Astronomical Union's [Minor Planet Center](https://www.minorplanetcenter.net).

based from a package developed by Chris Klein and Duncan Galloway