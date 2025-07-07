from resistify.coconat import coconat
from resistify.nlrexpress import nlrexpress
from resistify.hmmsearch import hmmsearch
from resistify.utility import parse_fasta, write_results
from resistify._loguru import logger


def nlr(args):
    sequences = parse_fasta(args.input)
    logger.info("Searching for NLRs...")
    sequences = hmmsearch(sequences, "nlr")
    if not args.retain:
        sequences = [sequence for sequence in sequences if sequence.has_nbarc]
        if len(sequences) == 0:
            logger.warning("No NLRs detected! Maybe try --retain?")
            return sequences
        else:
            logger.info(f"{len(sequences)} NLRs identified...")
    else:
        logger.info("NLRexpress will be run against all input sequences...")

    # If not specified, run NLRexpress in batches of 5 on all sequences or 1 on retained NLRs
    if args.retain and args.chunksize is None:
        chunksize = 5
    elif args.chunksize is None:
        chunksize = 1
    else:
        chunksize = args.chunksize

    sequences = nlrexpress(sequences, "all", chunksize, args.threads)

    if args.coconat:
        logger.info("Running CoCoNat to identify additional CC domains...")
        sequences = coconat(sequences)

    logger.info("Classifying NLRs...")

    for sequence in sequences:
        sequence.identify_lrr_domains(args.lrr_gap, args.lrr_length)
        if args.coconat:
            sequence.identify_cc_domains()
        sequence.merge_annotations(args.duplicate_gap)
        sequence.classify_nlr()

    write_results(sequences, args)
    return
