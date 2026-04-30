# Notebooks

## Asteroid Position Accuracy Notebook

Evaluates the accuracy of asteroid position predictions made by `pympc` by comparing them to online positions from SkyBoT, as a function of time separation from the catalogue epoch.

## NEA Position Accuracy Notebook

Quantifies how the accuracy of NEA (Near-Earth Asteroid) positions predicted by `pympc` changes as a function of time offset from the catalogue’s osculation epoch.

## Topocentric Corrections Notebook

Demonstrates the recovery of online topocentric positions from the Minor Planet Center by `pympc`'s own topocentric correction implementation, and quantifies the improvement in matching accuracy when applying topocentric corrections to geocentric `xephem` positions.
