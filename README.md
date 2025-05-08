<a href="https://ascl.net/2501.002"><img src="https://img.shields.io/badge/ascl-2501.002-blue.svg?colorB=262255" alt="ascl:2501.002" /></a>

pympc
=====

Perform checks for the presence of moving Solar-system bodies at astronomical locations for a given epoch.

## Installation

`pip install pympc`

or download/clone source and:

`python setup.py install`

## Setup

First, if we want to search for minor bodies, we need to grab the orbital elements catalogue from the Minor
Planet Center. This must be done at least once prior to any searches and can be run to overwrite the catalogues 
with the latest versions. The default call signature is shown.
```python
import pympc
xephem_cat =pympc.update_catalogue()
print(xephem_cat)
# e.g. /tmp/mpcorb_xephem.csv
```

The catalogue downloaded will be the [`mpcorb`](https://www.minorplanetcenter.net/data) catalogue.

The Near Earth Asteroid and Comets catalogues will be downloaded and used to update the `mpcorb` entries based on 
the values of the `include_nea` and `include_comets` arguments (both default to `True`).
 
It will create a csv file for each catalogue downloaded in the 
[xephem database format](http://www.clearskyinstitute.com/xephem/help/xephem.html#mozTocId468501) and return
the filepath to this file. By default, the file will be saved in the user's temporary directory - this can
be changed by setting the `cat_dir` argument.

Downloading this catalogue isn't necessary fi you just want to do major body checks (i.e. planets, moons.


## Usage 

Having downloaded the catalogue (see [Setup](#Setup)), we can now search for both major and minor bodies 
at a given location.

### Interactive searching

> **Note:** All information is output in logging. If you do not have a lgger set up in a session, running:
> ```python
> import logging
> fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
> logging.basicConfig(format=fmt, level=logging.INFO)
> ```
> prior to the examples will show this information.

Within an interpretor session, define a search location, epoch and radius and run the search.

```python
import astropy.units as u
import pympc
from astropy.time import Time

ra = 230.028 * u.deg
dec = -11.774 * u.deg
epoch = Time("2019-01-01T00:00")
search_radius = 5 * u.arcmin
observatory = 950  # observatory code for La Palma
pympc.minor_planet_check(ra, dec, epoch, search_radius, observatory=observatory)
```

Results are returned as an astropy table.

The above example uses explicit quantities, but if passed simple float arguments, and the program will assume the 
units (see comments below and `pympc.minor_planet_check()` docstring for unit assumptions).
```python
import pympc
ra = 230.028  # assumed degrees
dec = -11.774  # assumed degrees
epoch = 58484.  # assumed MJD
search_radius = 30  # assumed arcseconds
observatory = 950  # observatory code for La Palma
pympc.minor_planet_check(ra, dec, epoch, search_radius, observatory=observatory)
```

By default, the search will use a default filepath for the catalogue. if the file has been moved - or a 
custom `cat_dir` was passed to `pympc.update_catalogue()` - then the filepath can be specified.

```python
import pympc
pympc.minor_planet_check(ra=230.028, dec=-11.774, epoch=58484., search_radius=30, xephem_filepath='/path/to/mpcorb_xphem.csv')
```

The search will by default search both major and minor bodies. These can be toggled via the boolean arguments
`include_minor_bodies` and `include_major_bodies` arguments.

### Defining an observer

By default, if the `observatory` argument is not passed, the program will return geocentric coordinates. However, for
relatively nearby objects like minor bodies, there can be signicant parallax introduced by the location of an observer
on the Earth's surface. For this reason it is crucial to pass either an 
[observatory code](https://www.minorplanetcenter.net/iau/lists/ObsCodes.html) or a tuple containing the observatory
information. See the documentation for `pympc.minor_planet_check()` for more details.

### Console script searching

Installation of the package will create a `minor_planet_check` script, which can be accessed
from the command line. The options follow the same as the interactive searching, and results
are displayed as a table. For help on the command line use:
```bash
minor_planet_check --help
```

> **Note:** It is not currently possible to pass a custom set of observatory coordinates to the script - 
> an existing observatory code must be passed.



### Major bodies (planets, moons)

Major bodies (those with a `ephem.[body]` object in the `pyephem` package) can be included in the search by adding 
the `include_major_bodies` argument to `minor_planet_check()` when doing interactive searches, or 
`--match-to-major-bodies` when using the console script.

Additionally, `pympc` can perform a Hill sphere check for a position, to ensure it is not within the gravitational
sphere of influence of a planet. This is done by calling `planet_hill_sphere_check()` interactively, or adding
`--hill_sphere_check` to the console script.

#### Example of a minor Jovian moon

```bash
minor_planet_check 69.122371 21.11505 60695.428680 --match-to-major-bodies -r 5 --hill-sphere-check
```

This object is not in the major body catalogue of `pyephem`, but does provide a match to Jupiter's Hill sphere. Output:

```
Major and Minor Planet Check:
No major or minor bodies found.
Planet Hill Sphere Check:
  name          ra               dec             separation      mag 
------- ----------------- ------------------ ------------------ -----
Jupiter 69.83640596092219 21.603370530186574 2970.0900089441457 -2.46
```

## Speed and multiprocessing
The major body check takes much less than one second. The minor body check should take of order a few seconds, 
depending on multiprocessing capabilities.

The private function which actually performs the calculation is `_minor_planet_check()` (note leading underscore).
This can be called directly, to avoid the overhead associated with converting input arguments in `minor_planet_check()`,
if you provide them directly as required (see `_minor_planet_check()` docstring). Note that in this case a list of 
tuples is returned, rather than an astropy table.

By default, the program does calculate positions of bodies in the catalogue multiprocessed. To switch this off set
`chunk_size = 0`, i.e.:

```python
import pympc
pympc.minor_planet_check(ra=230.028, dec=-11.774, epoch=58484., search_radius=30, chunk_size=0)
```

## Limitations

1. The orbits are propagated following [xephem](http://www.clearskyinstitute.com/xephem) (via the 
[pyephem](https://rhodesmill.org/pyephem/) package), and this does not account for perturbations of the orbits. Thus, 
the accuracy of the position is dependent on the time difference between the epoch of the orbit elements and the epoch 
at which the search is being performed. Epoch differences between orbital elements calculation and observation of 
a few months or less will provide typical positional accuracies of less than a few arcsecond for the vast
majority of minor bodies. If the Minor Planet Center orbit epoch becomes many months out of date, then large separations
are to be expected. Note, additinoally, that a small number of bodies (those undergoing strong perturbations and
close to Earth) may be quite inaccurate (arcminutes) even at modest time differences between the search and orbit 
elements epochs. A fuller analysis is given in 
[notebooks/positional_accuracy.ipynb](notebooks/positional_accuracy.ipynb), with the following histogram showing the results.
![histogram showing positional accuracy of pympc vs minor planet center](/notebooks/position_accuracy.png "Histogram showing positional accuracy of `pympc` vs Minor Planet Center")
> **Note:** For this reason, `pympc` is not suitable to historical searches of positions since the Minor Planet Center 
> do not make available historical orbit elements catalgoues. It is intended for near real-time searches.

2. The `xephem` package, used to calculate positions, can only provide geocentric astrometric positions. `pympc` will 
calculate the topocentric correction as a post-processing to the initial position. The simple topometric correction 
applied is more than sufficient for the overwhelming majority of minor bodies, but for some near earth objects the 
correction can be large and the relatively simple treatment by `pympc` may not be sufficient. Additionally, in order 
to find matches in geocentric positions prior to applying the topocentric correction, a buffer is added to the search 
radius - this should capture the vast majority of cases where the geocentric position is outside the seach radius but 
the topocentric position is within it - unless the object is within ~1/3 AU of Earth. To work around this you can 
artifically inflate your search radius and filter yourself afterwards.

3. The filtering of matches based on magnitude via `max_mag` argument to `minor_planet_check()` is limited by the 
accuracy of the magnitude information in the database so some buffer should be applied to the desired magnitude cutoff 
to allow for this.

### Acknowledgements
This package makes use of data and/or services provided by the International Astronomical Union's 
[Minor Planet Center](https://www.minorplanetcenter.net).

Orbit elements are also sourced from [Lowell Observatory](https://asteroid.lowell.edu/main/), which is funded by the 
Lowell Observatory Endowment and NASA PDART grant NNX16AG52G.

Based from a package developed by Chris Klein and Duncan Galloway.
