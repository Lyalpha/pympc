.PHONY: help lint format type test all check-all build publish demo

help:
	@echo "Usage:"
	@echo "  make lint               Lint code with ruff (all files)"
	@echo "  make format             Format code with ruff (all files)"
	@echo "  make type               Type check with ty (all files)"
	@echo "  make test               Run unit tests"
	@echo "  make check-all          Run all checks (lint, format, type, test)"
	@echo "  make prek               Run prek (=pre-commit) hooks on all files"
	@echo "  make demo               Regenerate docs/cli_demo.svg for the README"
	@echo "  make build              Build the package wheel"
	@echo "  make publish            Publish package to PyPI"

lint:
	uv run ruff check pympc scripts test --fix

format:
	uv run ruff format pympc scripts test

type:
	uv run ty check pympc scripts test

test:
	uv run python -m unittest -v

check-all: lint format type test

prek:
	uv run prek run --all-files

demo:
	uv run python scripts/generate_cli_demo.py

build:
	uv build

publish: build
	uv publish
