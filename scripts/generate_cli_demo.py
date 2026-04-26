#!/usr/bin/env python
"""
Generate an SVG screenshot of the pympc CLI for use in README.md.

Run with:
    uv run python scripts/generate_cli_demo.py

The output SVG is written to docs/cli_demo.svg.
"""

import os
import subprocess
from pathlib import Path

from rich.console import Console
from rich.text import Text

REPO_ROOT = Path(__file__).parent.parent
OUT_DIR = REPO_ROOT / "docs" / "assets"
OUT_PATHS = [
    OUT_DIR / "cli_demo.svg",
    OUT_DIR / "jovian_moon_demo.svg",
]
COMMANDS = (
    [
        ["pympc", "--help"],
        ["pympc", "update-catalogue", "--help"],
        ["pympc", "update-obscode-cache", "--help"],
        ["pympc", "check", "--help"],
    ],
    [
        ["pympc", "check", "all", "69.122371", "21.11505", "60695.428680"],
    ],
)

WIDTH = 110


def _pympc_cmd(args: list[str]) -> str:
    """Invoke the pympc entry point directly, forcing ANSI colour output."""
    result = subprocess.run(
        ["uv", "run", "pympc"] + args,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        env={
            **os.environ,
            "FORCE_COLOR": "1",
            "COLUMNS": str(WIDTH),
        },
    )
    return result.stdout or result.stderr


def main() -> None:
    console = Console(record=True, width=WIDTH)

    for command, output_path in zip(COMMANDS, OUT_PATHS):
        for sub_command in command:
            prompt = Text(f"$ {' '.join(sub_command)}", style="bold green")
            console.print(prompt)
            output = _pympc_cmd(sub_command[1:])
            console.out(output, end="")
            console.print()

        output_path.parent.mkdir(parents=True, exist_ok=True)
        svg = console.export_svg(title="pympc CLI")
        output_path.write_text(svg)
        print(f"SVG written to {output_path}")


if __name__ == "__main__":
    main()
