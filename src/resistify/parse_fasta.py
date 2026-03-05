from pathlib import Path
import logging
from resistify.annotation import Protein

logger = logging.getLogger(__name__)


def parse_fasta(file_path: Path, sort: bool = True) -> dict[str, Protein]:
    logger.info(f"Loading sequences from {file_path}")

    def read_fasta(path):
        with open(path, "r") as f:
            header = None
            seq_parts = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        yield header, "".join(seq_parts)
                    header = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if header:
                yield header, "".join(seq_parts)

    proteins = {}

    for id, seq in read_fasta(str(file_path)):
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

    logger.info(f"{len(proteins)} sequences have been loaded")
    return proteins
