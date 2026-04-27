#!/usr/bin/env python
"""
Benchmark pympc.minor_planet_check with timing across chunk sizes.

Run with:
    uv run python scripts/benchmark.py
"""

import argparse
import time

import pympc

RA = 45.4625
DEC = 13.0417
EPOCH = "2026-04-26T17:18:00"
SEARCH_RADIUS = 50.0  # arcsec
OBSERVATORY = "La Palma"


def run_once(chunk_size: int) -> tuple[float, int]:
    from datetime import datetime

    epoch_dt = datetime.fromisoformat(EPOCH)
    t0 = time.perf_counter()
    result = pympc.minor_planet_check(
        ra=RA,
        dec=DEC,
        epoch=epoch_dt,
        search_radius=SEARCH_RADIUS,
        max_mag=None,
        include_minor_bodies=True,
        include_major_bodies=False,
        observatory=OBSERVATORY,
        chunk_size=chunk_size,
    )
    elapsed = time.perf_counter() - t0
    return elapsed, len(result)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark pympc.minor_planet_check")
    parser.add_argument(
        "--chunk-sizes",
        nargs="+",
        type=int,
        default=[0, 500, 1000, 5000, 20000, 100000],
        metavar="N",
        help="Chunk sizes to benchmark (0 = no multiprocessing)",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Number of repeats per chunk size (median is reported)",
    )
    args = parser.parse_args()

    print(
        f'Benchmark: RA={RA} Dec={DEC} epoch={EPOCH} radius={SEARCH_RADIUS}" obs={OBSERVATORY!r}'
    )
    print(f"Repeats per chunk size: {args.repeats}")
    print()
    print(
        f"{'chunk_size':>12}  {'min (s)':>9}  {'median (s)':>10}  {'max (s)':>9}  {'matches':>8}"
    )
    print("-" * 58)

    for chunk_size in args.chunk_sizes:
        label = str(chunk_size) if chunk_size > 0 else "0 (serial)"
        times = []
        nmatches = None
        for _ in range(args.repeats):
            elapsed, n = run_once(chunk_size)
            times.append(elapsed)
            nmatches = n
        times.sort()
        mid = len(times) // 2
        print(
            f"{label:>12}  {times[0]:>9.2f}  {times[mid]:>10.2f}  {times[-1]:>9.2f}  {nmatches:>8}"
        )

    print()


if __name__ == "__main__":
    main()
