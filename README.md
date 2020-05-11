pympc
=====

perform checks for the presence of minor bodies at astronomical locations for a given epoch.

### installation

`pip install pympc`

or download source and:

`python setup.py install`

### setup
frst import the package and grab the orbit element catalogue - this must be done at least once prior to doing any searches
and can be run to overwrite the catalogues with the latest versions. the default call signature is shown
```
>>> import pympc
>>> pympc.update_catalogue()
```

by default the catalogues downloaded will be the [`mpcorb`](https://www.minorplanetcenter.net/data) catalogue 
from the Minor Planet Center. [`astorb`](https://asteroid.lowell.edu/main/astorb) from Lowell Observatory may be
downloaded instead or as well. this behaviour can be changed using the `cat` argument.

the Near Earth Asteroid and Comets catalogues will be downloaded and used to update the `mpcorb` entries based on 
the values of the `include_nea` and `include_comets` arguments (both default to `True`).
 
 it will create a csv file for each catalogue downloaded in the 
 [xephem database format](http://www.clearskyinstitute.com/xephem/help/xephem.html#mozTocId468501).

##### example 1 (searching)
define our search location, epoch and radius and run the check
```
>>> import astropy.units as u
>>> from astropy.time import Time
>>> ra = 230.028 * u.deg
>>> dec = -11.774 * u.deg
>>> epoch = Time("2019-01-01T00:00")
>>> search_radius = 0.5 * u.arcmin
>>> pympc.minor_planet_check(ra, dec, epoch, search_radius)
```


##### example 2 (assumed units)
here we use float arguments, and the program assumes the units (see `pympc.minor_planet_check()` docstring)
```
>>> ra = 230.028  # assumed degrees
>>> dec = -11.774  # assumed degrees
>>> epoch = 58484.  # assumed MJD
>>> search_radius = 30  # assumed arcseconds
>>> pympc.minor_planet_check(ra, dec, epoch, search_radius)
```

##### example 3 (using a specific catalogue)
to use a specific orbit catalogue, specify the location of the xephem database , e.g.:

```
minor_planet_check(ra, dec, epoch, search_radius, xephem_filepath='/path/to/astorb_xphem.csv')
```

if omitted, it will attempt to use `/tmp/mpcorb_xephem.csv`.

### speed and multiprocessing
the check should take of order a second or two, depending on multiprocessing capabilities.

the private function which actually performs the calculation is `_minor_planet_check()` (note leading underscore).
this can be called to avoid the overhead associated with converting input arguments to `minor_planet_check()`, if
you can provide them directly as required (see `_minor_planet_check()` docstring).

by default the program will calculate positions of bodies in the catalogue multiprocessed. to switch this off set
`chunk_size = 0`, i.e.:

```
>>> pympc.minor_planet_check(ra, dec, epoch, search_radius, chunk_size=0)
```

### limitations
the orbits are propagated following [xephem](http://www.clearskyinstitute.com/xephem) (via the 
[pyephem](https://rhodesmill.org/pyephem/) package), and this does not account for perturbations of the orbits. thus 
the accuracy of the position is dependant on the time difference between the epoch of the orbit elements and the epoch 
at which the search is being performed. epoch differences between orbital elements calculation and observation of 
around a month or two should be fine for typical positional accuracies of a few arcsecond for most minor bodies - note
however that a small number of bodies (those under going strong perturbations) may be quite inaccurate (arcminutes).

the `xephem` package can only provide geocentric astrometric positions. as such, parallax effects for near-earth 
bodies will be significant, in addition to the lack of perturbation calculations above.

currently the epoch of the orbit elements is visible in the xephem db strings returned by `minor_planet_check()` as a
decimal year format (e.g. ..,2019.317808,..). some diagnostic information and warning when using large time differences
is to be implemented.

the filtering of matches based on magnitude via `max_mag` argument to `minor_planet_check()` is limited by the accuracy 
of the magnitude information in the database so some buffer should be applied to the desired magnitude cutoff to allow 
for this.

### acknowledgments
this package makes use of data and/or services provided by the International Astronomical Union's 
[Minor Planet Center](https://www.minorplanetcenter.net).

orbit elements are also sourced from [Lowell Observatory](https://asteroid.lowell.edu/main/), which is funded by the 
Lowell Observatory Endowment and NASA PDART grant NNX16AG52G.


based from a package developed by Chris Klein and Duncan Galloway.