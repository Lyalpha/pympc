# scripts/

Developer and maintenance scripts for `pympc`. None of these are part of the
installed package — they are run directly from the repo root, typically via 
`make` commands (see `Makefile`). Further details are provided below.

---

## `benchmark.py`

Benchmarks `pympc.minor_planet_check` wall-clock time across a range of
`chunk_size` values using a fixed sky position and epoch. Runs each
configuration a configurable number of times and reports min / median / max.

```bash
uv run python scripts/benchmark.py
uv run python scripts/benchmark.py --chunk-sizes 0 5000 20000 --repeats 5
# or simply
make benchmark
```

**Options**

| Flag | Default                             | Description |
|------|-------------------------------------|-------------|
| `--chunk-sizes N [N …]` | `0 500 1000 2000 10000 20000 50000` | Chunk sizes to test (`0` = serial) |
| `--repeats N` | `3`                                 | Runs per chunk size |

**Latest results (2026-04-27)**

```
Benchmark: RA=45.4625 Dec=13.0417 epoch=2026-04-26T17:18:00 radius=50.0" obs='La Palma'
Repeats per chunk size: 5

  chunk_size    min (s)  median (s)    max (s)   matches
----------------------------------------------------------
  0 (serial)      15.69       16.18      18.97         1
         500       9.91       12.31      13.54         1
        1000       9.93       10.43      11.01         1
        2000      10.21       11.20      11.92         1
       10000       9.12        9.79      12.20         1
       20000      10.15       11.12      11.54         1
       50000       8.68       10.24      11.76         1
```

-->

---

## `profile.py`

Minimal call to `minor_planet_check` (serial, `chunk_size=0`) so
it can be profiled without benchmark overhead. Intended for use with
`cProfile`.

```bash
# cProfile (two steps — or just run both via make)
uv run python -m cProfile -o scripts/prof.out scripts/profile.py
uv run python -c "import pstats; p=pstats.Stats('scripts/prof.out'); p.sort_stats('cumtime').print_stats(30)"
# or simply (runs both steps)
make profile
```

---

## `generate_cli_demos.py`

Renders SVG terminal captures of live `pympc` CLI output for use in
`docs/assets/` and the project README. Requires the catalogue to be up to
date (run `pympc update-catalogue` first).

Produces three SVGs:

| File | Content                                                                                 |
|------|-----------------------------------------------------------------------------------------|
| `docs/assets/cli_demo.svg` | All `--help` pages                                                                      |
| `docs/assets/jovian_moon_demo.svg` | Edge case check of a Jovian moon not in major bodies, but in the Hill sphere of Jupiter |
| `docs/assets/matches_demo.svg` | `check major` and `check minor` examples with real matches                              |

```bash
uv run python scripts/generate_cli_demos.py
# or
make demo
```

