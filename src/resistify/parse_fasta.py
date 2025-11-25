from pyfastx import Fastx
from pathlib import Path
import logging
from resistify.annotation import Protein

logger = logging.getLogger(__name__)


def parse_fasta(file_path: Path, sort: bool = True) -> dict[str, Protein]:
    logger.info(f"Loading sequencing from {file_path}")
    proteins = {}
    for id, seq in Fastx(file_path):
        seq = seq.strip("*")
        seq = seq.strip(".")

        if "*" in seq or "." in seq:
            logger.warning(
                f"{id} contains invalid * or . characters and will be skipped."
            )
            continue

        proteins[id] = Protein(id, seq)

    # Sort by largest first
    if sort:
        proteins = dict(
            sorted(
                proteins.items(), key=lambda item: len(item[1].sequence), reverse=True
            )
        )

    logger.info(f"{len(proteins)} sequences have been loaded from {file_path}")
    return proteins
