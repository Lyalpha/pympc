#!/usr/bin/env python
"""
Profiling for pympc.minor_planet_check.

Usage:
    uv run python -m cProfile -o scripts/prof.out scripts/profile.py
    uv run python -c "import pstats; p=pstats.Stats('scripts/prof.out'); p.sort_stats('cumtime').print_stats(30)"

"""

from datetime import datetime

import pympc

pympc.minor_planet_check(
    ra=45.4625,
    dec=13.0417,
    epoch=datetime.fromisoformat("2026-04-26T17:18:00"),
    search_radius=50.0,
    observatory="La Palma",
    include_minor_bodies=True,
    include_major_bodies=False,
    chunk_size=0,
)
