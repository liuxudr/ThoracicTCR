"""ThoracicTCR — pan-thoracic TCR repertoire toolkit command-line interface.

The CLI is organised as a tree of sub-apps so each pipeline stage has its own
namespace. Run ``thoracictcr --help`` to see all groups, or ``thoracictcr
GROUP --help`` to see the commands within one group.

Groups
------
* ``manifest``  Build TCGA / GEO sample manifests + atlas selection
* ``trust4``    TCR reconstruction (run-bam / run / fetch-refs / check)
* ``download``  Public TCR-antigen DB mirrors (VDJdb / McPAS / IEDB)
* ``metrics``   Diversity / clonality / VJ / clone-size / shared / rarefaction
* ``atlas``     Build and inspect the cross-cohort atlas joblib
* ``figures``   Render manuscript figures via Rscript
* ``info``      Environment + paths + tool version reporting
"""

from __future__ import annotations

import typer
from rich.console import Console

from . import __version__
from ._cli.atlas_app import app as atlas_app
from ._cli.download_app import app as download_app
from ._cli.figures_app import app as figures_app
from ._cli.info_app import app as info_app
from ._cli.manifest_app import app as manifest_app
from ._cli.metrics_app import app as metrics_app
from ._cli.trust4_app import app as trust4_app

app = typer.Typer(
    name="thoracictcr",
    help="Pan-thoracic TCR repertoire atlas + prognostic/ICI-response toolkit.",
    no_args_is_help=True,
    add_completion=False,
)
console = Console()

app.add_typer(manifest_app, name="manifest")
app.add_typer(trust4_app,   name="trust4")
app.add_typer(download_app, name="download")
app.add_typer(metrics_app,  name="metrics")
app.add_typer(atlas_app,    name="atlas")
app.add_typer(figures_app,  name="figures")
app.add_typer(info_app,     name="info")


@app.command("version")
def version() -> None:
    """Print the installed thoracictcr version."""
    console.print(f"thoracictcr [bold cyan]{__version__}[/bold cyan]")


if __name__ == "__main__":  # pragma: no cover
    app()
