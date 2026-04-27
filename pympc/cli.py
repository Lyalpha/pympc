from datetime import datetime
from numbers import Real

import rich_click as click


def _verbosity_to_level(verbose_count: int) -> str:
    return {0: "WARNING", 1: "INFO"}.get(min(verbose_count, 2), "DEBUG")


def _parse_epoch(epoch_value: str):
    try:
        return float(epoch_value)
    except ValueError:
        try:
            return datetime.fromisoformat(epoch_value)
        except ValueError as exc:
            raise click.BadParameter(
                "epoch must be either an MJD float (e.g. 60695.42868) or an ISO datetime "
                "(e.g. 2019-01-01T00:00:00)"
            ) from exc


def _effective_verbose(ctx, verbose_option: int) -> int:
    return max(int(ctx.obj.get("verbose", 0)), int(verbose_option))


def _parse_observatory(obs_string: str):
    """
    Parse observatory specification from CLI input.

    Accepts:
    - Integer obscode (e.g., "500")
    - 3-character string code (e.g., "950")
    - Observatory name (e.g., "La Palma")
    - Comma-separated tuple (e.g., "10.0,0.5,0.5")
    """
    obs_string = obs_string.strip()

    # Try comma-separated tuple first
    if "," in obs_string:
        parts = obs_string.split(",")
        if len(parts) != 3:
            raise click.BadParameter(
                f"Observatory tuple must have exactly 3 components (lon, rho_cos_phi, rho_sin_phi), got {len(parts)}"
            )
        try:
            return tuple(float(p.strip()) for p in parts)
        except ValueError as exc:
            raise click.BadParameter(
                f"Could not parse observatory tuple '{obs_string}' as floats: {exc}"
            ) from exc

    # Try integer conversion
    try:
        return int(obs_string)
    except ValueError:
        pass

    # Otherwise return as string (will be matched by code/name in get_observatory_data)
    return obs_string


def _render_astropy_table(title: str, table, hide_xephem: bool = True) -> None:
    from rich.console import Console
    from rich.table import Table

    console = Console()
    console.print(f"[bold]{title}[/bold]")

    if hide_xephem and "xephem_str" in table.colnames:
        table = table.copy()
        table.remove_column("xephem_str")

    if len(table) == 0:
        console.print("[dim]No matches found.[/dim]")
        return

    rich_table = Table(show_header=True, header_style="bold cyan")
    right_align = {"ra", "dec", "separation", "mag"}
    for col in table.colnames:
        justify = "right" if str(col) in right_align else "left"
        rich_table.add_column(str(col), justify=justify)

    def _fmt(col_name: str, value) -> str:
        if value is None:
            return ""
        if isinstance(value, Real):
            if col_name in {"ra", "dec"}:
                return f"{float(value):.6f}"
            if col_name == "separation":
                return f"{float(value):.3f}"
            if col_name == "mag":
                return f"{float(value):.2f}"
            return f"{float(value):.8g}"
        return str(value)

    for row in table:
        rich_table.add_row(*[_fmt(str(col), row[col]) for col in table.colnames])

    console.print(rich_table)


@click.group()
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity level of output. Add once for INFO, twice for DEBUG.",
)
@click.pass_context
def cli(ctx, verbose):
    """pympc command line interface."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


@cli.command("update-catalogue")
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity level of output. Add once for INFO, twice for DEBUG.",
)
@click.option(
    "--nea/--no-nea",
    "include_nea",
    default=True,
    show_default=True,
    help="Whether to include Near-Earth Asteroids (NEAs) in the catalogue.",
)
@click.option(
    "--comets/--no-comets",
    "include_comets",
    default=True,
    show_default=True,
    help="Whether to include comets in the catalogue.",
)
@click.option(
    "--cleanup/--no-cleanup",
    default=True,
    show_default=True,
    help="Whether to remove intermediate catalogue files before writing the final xephem database.",
)
@click.option(
    "--progress/--no-progress",
    "show_progress",
    default=True,
    show_default=True,
    help="Whether to show a progress bar during catalogue update.",
)
@click.option(
    "--catalogue-source",
    type=click.Choice(["astorb", "mpcorb"], case_sensitive=False),
    default="astorb",
    show_default=True,
    help="Source of the orbital catalogue data to use (AstOrb or MPCORB).",
)
@click.option(
    "--cat-dir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=str),
    help="Directory to store the final database. If not specified, the temporary directory will be used.",
)
@click.pass_context
def update_catalogue_cmd(
    ctx,
    verbose,
    include_nea,
    include_comets,
    cleanup,
    show_progress,
    catalogue_source,
    cat_dir,
):
    """Download/refresh orbital catalogues and generate an xephem CSV database."""
    from rich.console import Console

    from .pympc import update_catalogue
    from .utils import add_logging

    effective_verbose = _effective_verbose(ctx, verbose)
    add_logging(_verbosity_to_level(effective_verbose))

    filepath = update_catalogue(
        include_nea=include_nea,
        include_comets=include_comets,
        cat_dir=cat_dir,
        cleanup=cleanup,
        catalogue_source=catalogue_source,
        show_progress=show_progress,
    )
    Console().print(f"[green]Catalogue generated:[/green] {filepath}")


@cli.command("update-obscode-cache")
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity level of output. Add once for INFO, twice for DEBUG.",
)
@click.pass_context
def update_obscode_cache_cmd(ctx, verbose):
    """Refresh the locally cached MPC observatory-code database."""
    from rich.console import Console

    from .utils import add_logging, update_obscode_cache

    effective_verbose = _effective_verbose(ctx, verbose)
    add_logging(_verbosity_to_level(effective_verbose))

    update_obscode_cache()
    Console().print("[green]Observatory-code cache updated.[/green]")


@cli.command("check")
@click.argument(
    "mode",
    type=click.Choice(["all", "minor", "major", "hillsphere"], case_sensitive=False),
    default="all",
    required=False,
)
@click.argument("ra", type=float, help="Right Ascension in decimal degrees (e.g. 69.122371).")
@click.argument("dec", type=float, help="Declination in decimal degrees (e.g. 21.11505).")
@click.argument("epoch", type=str, help="Epoch of the position, either as MJD (e.g. 60695.42868) or ISO datetime (e.g. 2019-01-01T00:00:00).")
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity level of output. Add once for INFO, twice for DEBUG.",
)
@click.option(
    "-r",
    "--radius",
    type=float,
    default=5.0,
    show_default=True,
    help="Search radius in arcseconds for minor/major body checks.",
)
@click.option(
    "--max-mag",
    type=float,
    default=None,
    help="Only include matches with magnitude <= this value.",
)
@click.option(
    "--xephem-filepath",
    type=click.Path(exists=True, dir_okay=False, path_type=str),
    help="Explicit xephem catalogue to search. If omitted, latest source-matched file is used.",
)
@click.option("--chunk-size", type=int, default=20000, show_default=True,
              help="Number of catalogue entries to process in each multiprocessing chunk. "
                   "Set to 0 to disable multiprocessing.")
@click.option(
    "--catalogue-source",
    type=click.Choice(["astorb", "mpcorb"], case_sensitive=False),
    default="astorb",
    show_default=True,
    help="Base source used when auto-selecting latest xephem file.",
)
@click.option(
    "--cat-dir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=str),
    help="Directory to search for generated xephem catalogues.",
)
@click.option(
    "--observatory",
    type=str,
    default="500",
    show_default=True,
    help="MPC Observatory code, name, or custom coordinates. "
    "Code: integer or 3-char string (e.g. 500, 950). "
    "Name: e.g. 'La Palma', 'Greenwich'. "
    "Custom: comma-separated (longitude_deg, rho_cos_phi, rho_sin_phi) e.g. '10.0,0.5,0.5'.",
)
@click.pass_context
def check_cmd(
    ctx,
    mode,
    ra,
    dec,
    epoch,
    verbose,
    radius,
    max_mag,
    xephem_filepath,
    chunk_size,
    catalogue_source,
    cat_dir,
    observatory,
):
    """Run minor/major/hill-sphere checks at an input position and epoch."""
    from rich.console import Console

    from .pympc import minor_planet_check, planet_hill_sphere_check
    from .utils import add_logging

    effective_verbose = _effective_verbose(ctx, verbose)
    add_logging(_verbosity_to_level(effective_verbose))

    epoch_value = _parse_epoch(epoch)
    observatory_parsed = _parse_observatory(observatory)

    console = Console()
    console.print("\n[bold magenta]Search Parameters:[/bold magenta]")
    console.print(f"  Position (RA, Dec): ({ra:.6f}°, {dec:.6f}°)")
    console.print(f"  Search Radius: {radius:.2f} arcsec")
    console.print(f"  Epoch: {epoch}")
    if isinstance(observatory_parsed, tuple):
        console.print(
            f"  Observatory: Custom ({observatory_parsed[0]:.4f}, {observatory_parsed[1]:.6f}, {observatory_parsed[2]:.6f})"
        )
    else:
        console.print(f"  Observatory: [bold cyan]{observatory_parsed}[/bold cyan]")
    mode_display = (
        mode.title() if mode != "all" else "All (Minor & Major & Hill Sphere)"
    )
    if mode.lower() != "hillshpere":
        if xephem_filepath:
            console.print(f"  Catalogue: {xephem_filepath}")
        else:
            cat_dir_display = cat_dir or "system temp dir"
            console.print(
                f"  Catalogue: latest [bold cyan]{catalogue_source.upper()}[/bold cyan] xephem file in {cat_dir_display}"
            )
        console.print(f"  Search Mode: [bold cyan]{mode_display}[/bold cyan]")
        if max_mag is not None:
            console.print(f"  Max Magnitude: [bold cyan]{max_mag}[/bold cyan] mag")
    console.print()

    mode = mode.lower()

    if mode in ("all", "minor", "major"):
        include_minor = mode in ("all", "minor")
        include_major = mode in ("all", "major")
        table = minor_planet_check(
            ra=ra,
            dec=dec,
            epoch=epoch_value,
            search_radius=radius,
            xephem_filepath=xephem_filepath,
            max_mag=max_mag,
            include_minor_bodies=include_minor,
            include_major_bodies=include_major,
            observatory=observatory_parsed,
            chunk_size=chunk_size,
            cat_dir=cat_dir,
            catalogue_source=catalogue_source,
        )
        body_type = "Minor and Major" if mode == "all" else mode.capitalize()
        title = f"{body_type} Body Check Results (radius={radius} arcsec)"
        _render_astropy_table(title, table)

    if mode in ("all", "hillsphere"):
        hill_table = planet_hill_sphere_check(
            ra=ra,
            dec=dec,
            epoch=epoch_value,
            observatory=observatory_parsed,
        )
        _render_astropy_table("Hill Sphere Check Results", hill_table)
