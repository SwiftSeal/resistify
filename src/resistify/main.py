import typer
import os
import platform
import torch
from typing_extensions import Annotated
from pathlib import Path
from resistify.annotation import classify_nlrs, save_results
from resistify.hmmer import hmmsearch
from resistify.nlr_esm import predict_lrr
from resistify.console import console

DEFAULT_DEVICE = str(
    torch.device(
        "mps"
        if torch.backends.mps.is_built()
        else "cuda"
        if torch.cuda.is_available()
        else "cpu"
    )
)

app = typer.Typer()


@app.command()
def main(
    fasta: Path,
    output_dir: Annotated[
        Path, typer.Option(help="Output directory for results.")
    ] = Path("."),
    device: Annotated[
        str,
        typer.Option(help="Torch device to be used - can be 'cpu', 'cuda', or 'mps'."),
    ] = DEFAULT_DEVICE,
    quick: Annotated[
        bool, typer.Option(help="Enable to skip all ESM predictions.")
    ] = False,
    batch_size: Annotated[
        int,
        typer.Option(
            help="Batch size for ESM predictions - reduce this value for lower memory usage."
        ),
    ] = 8,
    pfam: Annotated[
        Path,
        typer.Option(
            help="Path to the Pfam database file - by default a smaller NLR-specific database is used."
        ),
    ] = Path(os.path.dirname(__file__)) / "data" / "db.hmm",
):
    console.rule("Resistify 2.0.0")
    console.log(f"Python version: {platform.python_version()}")
    console.log(f"PyTorch version: {torch.__version__}")
    console.log(f"Using device: {device}")
    console.log(f"OS: {platform.system()} {platform.machine()}")

    proteins = hmmsearch(fasta, pfam)

    if not quick:
        predict_lrr(proteins, device, batch_size=batch_size)

    classify_nlrs(proteins)

    save_results(proteins, output_dir)

    console.print("Thank you for using Resistify!", style="bold")
    console.print(":page_facing_up: https://doi.org/10.1177/11779322241308944")
