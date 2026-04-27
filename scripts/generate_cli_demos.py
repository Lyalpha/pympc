#!/usr/bin/env python
"""
Generate SVG terminal captures of pympc CLI output for README/docs.

Run with:
    uv run python scripts/generate_cli_demos.py
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
    OUT_DIR / "matches_demo.svg",
]
COMMANDS = (
    [
        ["pympc", "--help"],
    ],
    [
        ["pympc", "check", "all", "69.122371", "21.11505", "60695.428680"],
    ],
    [
        ["pympc", "check", "major", "8.442083", "1.290639", "2026-04-26T12:00:00"],
        ["pympc", "check", "minor", "214.768413", "10.41659", "61157.006230"],
        ["pympc", "check", "minor", "84.2229166", "21.5972222", "61157.66100"],
    ],
)

WIDTH = 110


def _pympc_cmd(args: list[str]) -> str:
    """Invoke the pympc entry point directly and return combined CLI output."""
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


def _render_block(commands: list[list[str]], output_path: Path) -> None:
    console = Console(record=True, width=WIDTH)

    for sub_command in commands:
        console.print(Text(f"$ {' '.join(sub_command)}", style="bold green"))
        output = _pympc_cmd(sub_command[1:])
        console.print(Text.from_ansi(output), end="")
        console.print()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(console.export_svg(title="pympc"))
    print(f"SVG written to {output_path}")


def main() -> None:
    for command_group, output_path in zip(COMMANDS, OUT_PATHS):
        _render_block(command_group, output_path)


if __name__ == "__main__":
    main()
