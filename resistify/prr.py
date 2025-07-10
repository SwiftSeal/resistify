from pathlib import Path
from resistify._loguru import logger


def prr(
    infile: Path,
    outdir: Path = Path("."),
    lrr_gap: int = 75,
    lrr_length: int = 4,
    duplicate_gap: int = 100,
    chunksize: int = 5,
    threads: int = 1,
):
    """
    Identify and classify PRR sequences.
    """
    try:
        from resistify.tmbed import tmbed
    except ImportError:
        raise ImportError(
            "Some TMbed dependencies are missing - install resistify[full] to use this feature."
        )
    from resistify.hmmsearch import hmmsearch
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
