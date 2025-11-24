import sys
import os
import logging
from pathlib import Path
import pyhmmer
from resistify.annotation import Protein, Annotation

logger = logging.getLogger(__name__)

NLR_HMM_DB = Path(os.path.dirname(__file__)) / "data" / "nlr.hmm"
RLP_HMM_DB = Path(os.path.dirname(__file__)) / "data" / "rlp.hmm"

NLR_ACCESSION_DOMAINS = {
    "PF05659": "RPW8",
    "PF18052": "CC",
    "PF01582": "TIR",
    "PF13676": "TIR",
    "PF00931": "NBARC",
    "PF00560": "LRR",
    "PF07725": "LRR",
    "PF12061": "LRR",
    "PF12799": "LRR",
    "PF13855": "LRR",
    "PF23211": "LRR",
    "PF23247": "LRR",
    "PF23286": "LRR",
    "PF23598": "LRR",
    "PF23622": "LRR",
    "PF23952": "LRR",
    "PF24758": "LRR",
    "PF25013": "LRR",
    "PF25019": "LRR",
    "PF20160": "C-JID",
}

RLP_ACCESSION_DOMAINS = {
    "PF00069": "PKinase",
    "PF01453": "G-LecRLK",
    "PF00954": "G-LecRLK",
    "PF08276": "G-LecRLK",
    "PF00024": "G-LecRLK",
    "PF14295": "G-LecRLK",
    "PF00139": "L-LecRLK",
    "PF00059": "C-LecRLK",
    "PF13947": "WAK",
    "PF14380": "WAK",
    "PF00008": "WAK",
    "PF08488": "WAK",
    "PF07645": "WAK",
    "PF12662": "WAK",
    "PF12947": "WAK",
    "PF11721": "CrRLK1L",
    "PF12819": "CrRLK1L",
    "PF01476": "LysM",
    "PF01657": "CRK",
    "PF00314": "Thaumatin",
    "PF13540": "CR-like",
    "PF19160": "SPARK",
    "PF00704": "GH18",
    "PF00182": "GH19",
    "PF00188": "CAP",
    "PF16101": "PRIMA1",
}

ACCESSION_DOMAINS = {**NLR_ACCESSION_DOMAINS, **RLP_ACCESSION_DOMAINS}


def hmmsearch(fasta: Path, database: Path, threads: int = 0) -> dict[str, Protein]:
    alphabet = pyhmmer.easel.Alphabet.amino()
    proteins = {}
    sequences = []

    with pyhmmer.easel.SequenceFile(fasta, digital=True, alphabet=alphabet) as seq_file:
        sequences = list(seq_file)
        for sequence in sequences:
            id = sequence.name.decode()
            aa = alphabet.decode(sequence.sequence)
            proteins[id] = Protein(id, aa)

    logger.info(f"{len(proteins)} sequences have been loaded")

    with pyhmmer.plan7.HMMFile(database) as hmm_file:
        for tophit in pyhmmer.hmmsearch(
            hmm_file, sequences, bit_cutoffs="gathering", Z=45638612, cpus=threads
        ):
            accession = tophit.query.accession.decode().split(".")[0]
            try:
                name = ACCESSION_DOMAINS[accession]
            except KeyError:
                name = tophit.query.name.decode()
            for hit in tophit:
                id = hit.name.decode()
                # Can use included to only add high-scoring domain matches
                for domain in hit.domains.included:
                    proteins[id].add_annotation(
                        Annotation(
                            name,
                            domain.env_from,
                            domain.env_to,
                            type="domain",
                            score=float(domain.score),
                            source="hmmer",
                            accession=accession,
                        )
                    )
    return proteins


def update_pfam_db(pfam_path: Path):
    with (
        pyhmmer.plan7.HMMFile(pfam_path) as hmm_file,
        open(NLR_HMM_DB, "wb") as nlr,
        open(RLP_HMM_DB, "wb") as rlp,
    ):
        for hmm in hmm_file:
            accession = hmm.accession.decode().split(".")[0]
            if accession in NLR_ACCESSION_DOMAINS.keys():
                hmm.write(nlr, binary=False)
            elif accession in RLP_ACCESSION_DOMAINS.keys():
                hmm.write(rlp, binary=False)


if __name__ == "__main__":
    update_pfam_db(Path(sys.argv[1]))
