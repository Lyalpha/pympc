.PHONY: help lint format type test all check-all build publish demo benchmark profile

help:
	@echo "Usage:"
	@echo "  make lint               Lint code with ruff (all files)"
	@echo "  make format             Format code with ruff (all files)"
	@echo "  make type               Type check with ty (all files)"
	@echo "  make test               Run unit tests"
	@echo "  make check-all          Run all checks (lint, format, type, test)"
	@echo "  make prek               Run prek (=pre-commit) hooks on all files"
	@echo "  make demo               Regenerate docs/assets SVGs for the README"
	@echo "  make benchmark          Benchmark minor_planet_check across chunk sizes"
	@echo "  make profile            Profile minor_planet_check with cProfile"
	@echo "  make build              Build the package wheel"
	@echo "  make publish            Publish package to PyPI"

lint:
	uv run ruff check pympc scripts test --fix

format:
	uv run ruff format pympc scripts test

type:
	uv run ty check pympc scripts test

test:
	uv run python -m pytest -v

check-all: lint format type test

prek:
	uv run prek run --all-files

demo:
	uv run python scripts/generate_cli_demos.py

benchmark:
	uv run python scripts/benchmark.py

profile:
	uv run python -m cProfile -o scripts/prof.out scripts/profile.py
	uv run python -c "import pstats; p=pstats.Stats('scripts/prof.out'); p.sort_stats('cumtime').print_stats(30)"

build:
	uv build

publish: build
	uv publish
