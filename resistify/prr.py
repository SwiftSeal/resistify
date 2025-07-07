from resistify.hmmsearch import hmmsearch
from resistify.tmbed import tmbed
from resistify.nlrexpress import nlrexpress
from resistify.utility import parse_fasta, write_results
from resistify._loguru import logger


def prr(args):
    if args.chunksize is None:
        chunksize = 5
    else:
        chunksize = args.chunksize

    logger.info("Searching for PRRs...")
    sequences = parse_fasta(args.input)
    sequences = hmmsearch(sequences, "prr")

    sequences = tmbed(sequences)
    sequences = [sequence for sequence in sequences if sequence.is_rlp()]
    if len(sequences) > 0:
        logger.info(f"{len(sequences)} PRRs identified...")
        sequences = nlrexpress(sequences, "lrr", chunksize, args.threads)

        logger.info("Classifying PRRs...")
        for sequence in sequences:
            sequence.identify_lrr_domains(args.lrr_gap, args.lrr_length)
            sequence.merge_annotations(args.duplicate_gap)
            sequence.classify_rlp()
    else:
        logger.warning("No PRRs detected!")

    write_results(sequences, args)
    return
