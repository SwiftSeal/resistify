import typer
import os
from pathlib import Path
import pyhmmer
from resistify.annotation import Protein, Annotation
from resistify.console import console

app = typer.Typer()

ACCESSION_DOMAINS = {
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
}


def hmmsearch(fasta: Path, pfam_path: Path) -> list[Protein]:
    with console.status("Running HMMER search..."):
        alphabet = pyhmmer.easel.Alphabet.amino()

        # Initialise a dict of protein sequences - convert to list at end
        proteins = {}
        with pyhmmer.easel.SequenceFile(fasta, digital=True, alphabet=alphabet) as seq_file:
            sequences = list(seq_file)
            for sequence in sequences:
                id = sequence.name.decode()
                aa = alphabet.decode(sequence.sequence)
                proteins[id] = Protein(id, aa)

        console.log(f"{len(proteins)} sequences have been loaded")

        with pyhmmer.plan7.HMMFile(pfam_path) as hmm_file:
            for tophit in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="gathering", Z=45638612):
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
                                score=float(domain.score),
                                accession=accession,
                            )
                        )
    return list(proteins.values())


@app.command()
def update_pfam_db(pfam_path: Path):
    """
    An internal script designed to use the core ACCESSION_DOMAINS to produce a small HMM file to be provided with code
    """
    outfile = Path(os.path.dirname(__file__)) / "data" / "db.hmm"
    with pyhmmer.plan7.HMMFile(pfam_path) as hmm_file, open(outfile, "wb") as f:
        for hmm in hmm_file:
            accession = hmm.accession.decode().split(".")[0]
            if accession in ACCESSION_DOMAINS.keys():
                hmm.write(f, binary=False)

if __name__ == "__main__":
    app()
