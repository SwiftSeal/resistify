import os
from pathlib import Path
import typer
from typing_extensions import Annotated
from resistify._loguru import logger
from resistify.utility import hello, goodbye, get_threads

app = typer.Typer(no_args_is_help=True)


@app.callback()
def callback(verbose: bool = False):
    if verbose:
        logger.update_level("DEBUG")

    hello()


@app.command(short_help="Identify and classify NLR sequences.")
def nlr(
    infile: Annotated[
        Path,
        typer.Argument(
            help="Input FASTA file containing sequences to analyze.",
            exists=True,
            readable=True,
        ),
    ],
    outdir: Annotated[
        Path,
        typer.Option(
            "--outdir", "-o", help="Directory to save results.", writable=True
        ),
    ] = Path("."),
    run_coconat: Annotated[
        bool,
        typer.Option(
            "--coconat", "-c", help="Run CoCoNat to identify additional CC domains."
        ),
    ] = False,
    retain: Annotated[
        bool,
        typer.Option(
            "--retain",
            "-r",
            help="Retain all sequences, not just those with NB-ARC domains.",
        ),
    ] = False,
    lrr_gap: Annotated[
        int,
        typer.Option(
            "--lrr-gap", help="Minimum gap (in amino acids) between LRR motifs."
        ),
    ] = 75,
    lrr_length: Annotated[
        int,
        typer.Option(
            "--lrr-length",
            help="Minimum number of LRR motifs required to be considered an LRR domain.",
        ),
    ] = 4,
    duplicate_gap: Annotated[
        int,
        typer.Option(
            "--duplicate-gap",
            help="Gap size (in amino acids) to consider merging duplicate annotations.",
        ),
    ] = 100,
    chunksize: Annotated[
        int | None,
        typer.Option(
            "--chunksize", help="Number of sequences per split for jackhmmer."
        ),
    ] = None,
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            "-t",
            help="Number of threads available for NLRexpress. Defaults to all available threads",
        ),
    ] = get_threads(),
):
    """
    Identify and classify NLR sequences.
    This command will search for NLRs using HMMs, run NLRexpress to identify NLR-associated motifs, and optionally run CoCoNat to identify additional CC domains.
    Together, this information will be used to classify NLRs into their respective classes.
    """
    from resistify.coconat import coconat
    from resistify.nlrexpress import nlrexpress
    from resistify.hmmsearch import hmmsearch
    from resistify.utility import parse_fasta, write_results

    sequences = parse_fasta(infile)
    logger.info("Searching for NLRs...")
    sequences = hmmsearch(sequences, "nlr")
    if not retain:
        sequences = [sequence for sequence in sequences if sequence.has_nbarc]
        if len(sequences) == 0:
            logger.warning("No NLRs detected! Maybe try --retain?")
            return sequences
        else:
            logger.info(f"{len(sequences)} NLRs identified...")
    else:
        logger.info("NLRexpress will be run against all input sequences...")

    # If not specified, run NLRexpress in batches of 5 on all sequences or 1 on retained NLRs
    if retain and chunksize is None:
        chunksize = 5
    elif chunksize is None:
        chunksize = 1
    else:
        chunksize = chunksize

    sequences = nlrexpress(sequences, "all", chunksize, threads)

    if run_coconat:
        logger.info("Running CoCoNat to identify additional CC domains...")
        sequences = coconat(sequences)

    logger.info("Classifying NLRs...")

    for sequence in sequences:
        sequence.identify_lrr_domains(lrr_gap, lrr_length)
        if run_coconat:
            sequence.identify_cc_domains()
        sequence.merge_annotations(duplicate_gap)
        sequence.classify_nlr()

    write_results(sequences, outdir, "nlr", coconat=run_coconat, retain=retain)
    goodbye(coconat=run_coconat)


@app.command()
def prr(
    infile: Annotated[
        Path,
        typer.Argument(
            help="Input FASTA file containing sequences to analyze.", exists=True
        ),
    ],
    outdir: Annotated[
        Path, typer.Option("--outdir", "-o", help="Directory to save results.")
    ] = Path("."),
    lrr_gap: Annotated[
        int,
        typer.Option(
            "--lrr-gap", help="Minimum gap (in amino acids) between LRR motifs."
        ),
    ] = 75,
    lrr_length: Annotated[
        int,
        typer.Option(
            "--lrr-length",
            help="Minimum number of LRR motifs required to be considered an LRR domain.",
        ),
    ] = 4,
    duplicate_gap: Annotated[
        int,
        typer.Option(
            "--duplicate-gap",
            help="Gap size (in amino acids) to consider merging duplicate annotations.",
        ),
    ] = 100,
    chunksize: Annotated[
        int,
        typer.Option(
            "--chunksize", help="Number of sequences per split for jackhmmer."
        ),
    ] = 5,
    threads: Annotated[
        int,
        typer.Option(
            "--threads", "-t", help="Number of threads available for NLRexpress."
        ),
    ] = get_threads(),
):
    """
    Identify and classify PRR sequences.
    """
    from resistify.hmmsearch import hmmsearch
    from resistify.tmbed import tmbed
    from resistify.nlrexpress import nlrexpress
    from resistify.utility import parse_fasta, write_results

    logger.info("Searching for PRRs...")
    sequences = parse_fasta(infile)
    sequences = hmmsearch(sequences, "prr")

    sequences = tmbed(sequences)
    sequences = [sequence for sequence in sequences if sequence.is_rlp()]
    if len(sequences) > 0:
        logger.info(f"{len(sequences)} PRRs identified...")
        sequences = nlrexpress(sequences, "lrr", chunksize, threads)

        logger.info("Classifying PRRs...")
        for sequence in sequences:
            sequence.identify_lrr_domains(lrr_gap, lrr_length)
            sequence.merge_annotations(duplicate_gap)
            sequence.classify_rlp()
    else:
        logger.warning("No PRRs detected!")

    write_results(sequences, outdir, "prr")
    goodbye(tmbed=True)


@app.command(short_help="Draw NLRs or PRRs from results directory.")
def draw(
    results_dir: Annotated[
        str, typer.Argument(help="Directory containing results files.")
    ],
    outfile: Annotated[
        Path, typer.Option("--outfile", "-o", help="Output file for the plot.")
    ] = Path("resistify_plot.png"),
    query: Annotated[
        str | None,
        typer.Option(
            "--query",
            "-q",
            help="Comma-separated list of sequence names to plot. If not provided, all sequences will be plotted.",
        ),
    ] = None,
    width: Annotated[
        int | None,
        typer.Option("--width", "-w", help="Width of the output figure in inches."),
    ] = 10,
    height: Annotated[
        int,
        typer.Option("--height", "-h", help="Height of the output figure in inches."),
    ] = 5,
    hide_motifs: Annotated[
        bool,
        typer.Option("--hide-motifs", help="Hide motifs in the plot.", is_flag=True),
    ] = False,
):
    """
    Draw NLR or PRR sequences from the results directory.
    This command will automatically detect whether to draw NLRs or PRRs based on the presence of `nlr.fasta` or `prr.fasta` in the results directory.
    I'd recommend using --query to specify which sequences to plot, as otherwise it will plot all sequences in the directory!
    """
    from resistify.draw import collect_data, draw_nlr, draw_prr
    import os

    if query is None:
        logger.info("No specific queries provided, plotting all sequences.")
        queries = None
    else:
        queries = [q.strip() for q in query.split(",") if q.strip()]

    sequence_data = collect_data(queries, results_dir)

    # First, need to detect if we are drawing NLRs or PRRs
    # Can do this by checking for prr.fasta or nlr.fasta in directory
    if "nlr.fasta" in os.listdir(results_dir):
        logger.info("Detected NLR sequences, drawing NLRs...")
        draw_nlr(sequence_data, outfile, width, height, hide_motifs)
    elif "prr.fasta" in os.listdir(results_dir):
        logger.info("Detected PRR sequences, drawing PRRs...")
        draw_prr(sequence_data, outfile, width, height, hide_motifs)


@app.command(short_help="Manually download ESM and ProtT5 models.")
def download_models():
    """
    Manually download models - by default, models will be downloaded automatically when required.
    This command is only necessary if you want to download the models in advance or if you are running in an environment without internet access.
    """
    from transformers import T5EncoderModel, T5Tokenizer
    import esm

    hf_home = os.environ.get("HF_HOME", "~/.cache/huggingface/")
    torch_home = os.environ.get("TORCH_HOME", "~/.cache/torch/")

    logger.info(f"Hugging Face models will be stored under {hf_home}")
    logger.info(f"Torch models will be stored under {torch_home}")

    # ProtT5
    logger.info("Loading ProtT5 models...")
    T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
    T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")

    # ESM
    logger.info("Loading ESM models...")
    esm.pretrained.esm2_t33_650M_UR50D()

    logger.info("Models have been downloaded successfully")


if __name__ == "__main__":
    app()
