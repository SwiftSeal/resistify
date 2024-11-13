import sys
import os
import csv
import logging
from Bio import SeqIO
from resistify.annotations import Sequence

log = logging.getLogger(__name__)


def create_output_directory(outdir):
    try:
        expanded_outdir = os.path.expanduser(os.path.expandvars(outdir))
        os.makedirs(expanded_outdir, exist_ok=True)
        log.debug(f"Output directory created at {expanded_outdir}")
        return expanded_outdir
    except OSError as e:
        log.error(f"Error creating output directory: {e}")
        sys.exit(1)


def parse_fasta(path):
    sequences = []
    with open(path) as file:
        for record in SeqIO.parse(file, "fasta"):
            # need to remove asterisk, interferes with hmmsearch
            sequence_str = str(record.seq).strip("*")
            if "*" in sequence_str:
                log.error(f"Internal stop codon detected in sequence {record.id}")
                sys.exit(1)
            if len(sequence_str) > 100000:
                log.error(f"Sequence {record.id} is too long (>100000 aa)")
                sys.exit(1)
            sequences.append(Sequence(record.id, sequence_str))

    return sequences


def save_fasta(sequences, path, nlr_only=False):
    with open(path, "w") as file:
        for sequence in sequences:
            if nlr_only and sequence.classification is None:
                continue
            file.write(f">{sequence.id}\n")
            file.write(f"{sequence.seq}\n")
    return path


def result_table(sequences, results_dir):
    with open(os.path.join(results_dir, "results.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(
            [
                "Sequence",
                "Length",
                "Motifs",
                "Domains",
                "Classification",
                "NBARC_motifs",
                "MADA",
                "MADAL",
                "CJID",
            ]
        )

        nbarc_motifs = [
            "VG",
            "P-loop",
            "RNSB-A",
            "Walker-B",
            "RNSB-B",
            "RNSB-C",
            "RNSB-D",
            "GLPL",
            "MHD",
        ]

        for sequence in sequences:
            n_nbarc_motifs = 0
            for motif in nbarc_motifs:
                if len(sequence.motifs[motif]) > 0:
                    n_nbarc_motifs += 1

            table_writer.writerow(
                [
                    sequence.id,
                    len(sequence.seq),
                    sequence.motif_string(),
                    sequence.domain_string,
                    sequence.classification,
                    n_nbarc_motifs,
                    sequence.mada,
                    sequence.madal,
                    sequence.cjid,
                ]
            )


def domain_table(sequences, results_dir):
    with open(os.path.join(results_dir, "domains.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(["Sequence", "Domain", "Start", "End"])
        for sequence in sequences:
            for annotation in sequence.merged_annotations:
                table_writer.writerow(
                    [
                        sequence.id,
                        annotation.domain,
                        annotation.start,
                        annotation.end,
                    ]
                )


def annotation_table(sequences, results_dir):
    with open(os.path.join(results_dir, "annotations.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(
            ["Sequence", "Domain", "Start", "End", "E_value", "Score", "Source"]
        )
        for sequence in sequences:
            for annotation in sequence.annotations:
                table_writer.writerow(
                    [
                        sequence.id,
                        annotation.domain,
                        annotation.start,
                        annotation.end,
                        annotation.evalue,
                        annotation.score,
                        annotation.source,
                    ]
                )


def motif_table(sequences, results_dir):
    from resistify.nlrexpress import MOTIF_SPAN_LENGTHS

    output_path = os.path.join(results_dir, "motifs.tsv")
    with open(output_path, "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(
            [
                "Sequence",
                "Motif",
                "Position",
                "Probability",
                "Downstream_sequence",
                "Motif_sequence",
                "Upstream_sequence",
            ]
        )
        for sequence in sequences:
            for motif in sequence.motifs:
                for item in sequence.motifs[motif]:
                    aa_sequence = sequence.seq
                    downstream_sequence = aa_sequence[item.position - 5 : item.position]
                    motif_sequence = aa_sequence[
                        item.position : item.position + MOTIF_SPAN_LENGTHS[motif]
                    ]
                    upstream_sequence = aa_sequence[
                        item.position
                        + MOTIF_SPAN_LENGTHS[motif] : item.position
                        + MOTIF_SPAN_LENGTHS[motif]
                        + 5
                    ]
                    table_writer.writerow(
                        [
                            sequence.id,
                            motif,
                            item.position,
                            item.probability,
                            downstream_sequence,
                            motif_sequence,
                            upstream_sequence,
                        ]
                    )


def extract_nbarc(sequences, results_dir):
    """
    Extract all nbarc domains of all proteins into a fasta file.
    """
    with open(os.path.join(results_dir, "nbarc.fasta"), "w") as file:
        for sequence in sequences:
            count = 1
            for annotation in sequence.merged_annotations:
                if annotation.domain == "NB-ARC":
                    file.write(f">{sequence.id}_{count}\n")
                    file.write(
                        f"{sequence.seq[annotation.start:annotation.end]}\n"
                    )
                    count += 1


def coconat_table(sequences, results_dir):
    output_path = os.path.join(results_dir, "coconat.tsv")
    with open(output_path, "w") as f:
        for sequence in sequences:
            for i, probability in enumerate(sequence.cc_probs):
                f.write(f"{sequence.id}\t{i}\t{probability}\n")
