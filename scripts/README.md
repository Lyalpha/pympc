# scripts/

Developer and maintenance scripts for `pympc`. None of these are part of the
installed package — they are run directly from the repo root, typically via
`make` commands (see `Makefile`).

They typically require the `pympc catalogue update` has been run prior to
invoking them, so that the xephem database has been created.

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
uv run python scripts/benchmark.py
Benchmark: RA=45.4625 Dec=13.0417 epoch=2026-04-26T17:18:00 radius=50.0" obs='La Palma'
Repeats per chunk size: 3

  chunk_size    min (s)  median (s)    max (s)   matches
----------------------------------------------------------
  0 (serial)      12.83       12.84      14.51         1
         500       5.39        5.61       7.22         1
        1000       4.82        5.08       7.86         1
        5000       3.92        3.98       4.24         1
       20000       3.93        3.96       4.07         1
      100000       4.28        4.40       4.47         1
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
