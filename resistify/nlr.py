from pathlib import Path
from resistify._loguru import logger


def nlr(
    infile: Path,
    outdir: Path,
    run_coconat: bool = False,
    retain: bool = False,
    lrr_gap: int = 75,
    lrr_length: int = 4,
    duplicate_gap: int = 100,
    chunksize: int | None = None,
    threads: int = 1,
):
    """
    Identify and classify NLR sequences.
    This command will search for NLRs using HMMs, run NLRexpress to identify NLR-associated motifs, and optionally run CoCoNat to identify additional CC domains.
    Together, this information will be used to classify NLRs into their respective classes.
    """
    if run_coconat:
        try:
            from resistify.coconat import coconat
        except ImportError:
            raise ImportError(
                "Some CoCoNat dependencies are missing - install resistify[full] to use this feature."
            )
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
