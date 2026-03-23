import csv
import logging
import re
import shutil
import tarfile
from pathlib import Path

import svg
from svg import Element

from resistify.annotation import Protein

logger = logging.getLogger(__name__)

COLOUR_PALETTE = {
    "CC": "#648FFF",
    "RPW8": "#785EF0",
    "TIR": "#DC267F",
    "NBARC": "#FE6100",
    "LRR": "#FFB000",
    "signal_peptide": "#648FFF",
    "alpha_inwards": "#DC267F",
    "PKinase": "#FE6100",
    "extEDVID": "#648FFF",
    "bA": "#DC267F",
    "aA": "#DC267F",
    "bC": "#DC267F",
    "aC": "#DC267F",
    "bDaD1": "#DC267F",
    "aD3": "#DC267F",
    "VG": "#FE6100",
    "P-loop": "#FE6100",
    "RNBS-A": "#FE6100",
    "Walker-B": "#FE6100",
    "RNBS-B": "#FE6100",
    "RNBS-C": "#FE6100",
    "RNBS-D": "#FE6100",
    "GLPL": "#FE6100",
    "MHD": "#FE6100",
    "LxxLxL": "#FFB000",
}

DISPLAY_NAMES = {
    "signal_peptide": "SP",
    "alpha_inwards": "TM",
    "alpha_outwards": "TM",
}


def draw_svg(protein: Protein, output_path: Path, command: str):
    elements: list[Element]

    elements = [
        svg.Line(
            x1=0,
            x2=protein.length,
            y1=60,
            y2=60,
            stroke="black",
            stroke_width=6,
        ),
    ]

    motifs = [a for a in protein.annotations if a.type == "motif"]

    for motif in motifs:
        elements.append(
            svg.Line(
                x1=motif.start,
                x2=motif.start,
                y1=32,
                y2=60,
                stroke="black",
                stroke_width=2,
            ),
        )
        elements.append(
            svg.Text(
                x=motif.start - 8,  # gentle shift
                y=28,
                text=motif.name,
                transform=[svg.Rotate(-45, motif.start, 20)],
                font_size=8,
                font_family="sans-serif",
                fill=COLOUR_PALETTE.get(motif.name, "black"),
            )
        )

    domains = [
        a for a in protein.annotations if a.type == "domain" and a.source == "merged"
    ]
    if command == "nlr":
        domains = [d for d in domains if d.name != "MADA"]
    if command == "prr":
        domains = [d for d in domains if d.name]
    for domain in domains:
        elements.append(
            svg.Rect(
                x=domain.start,
                y=45,
                width=domain.end - domain.start + 1,
                height=30,
                fill=COLOUR_PALETTE.get(domain.name, "grey"),
                stroke="black",
                stroke_width=2,
            )
        )
        elements.append(
            svg.Text(
                x=domain.start + (domain.end - domain.start) / 2,
                y=60,
                text=DISPLAY_NAMES.get(domain.name, domain.name),
                font_size=14,
                font_family="sans-serif",
                text_anchor="middle",
                dominant_baseline="middle",
            )
        )

    canvas = svg.SVG(width=protein.length, height=100, elements=elements)

    with open(output_path, "w") as f:
        f.write(str(canvas))


def save_results(
    proteins: dict[str, Protein], output_dir: Path, command: str, draw: bool = True
):
    logger.info("Saving results")
    output_dir.mkdir(parents=True, exist_ok=True)
    if draw:
        (output_dir / "plots").mkdir(parents=True, exist_ok=True)

    classified_fasta_path = output_dir / f"{command}.fa"

    with (
        open(output_dir / "results.tsv", "w", newline="") as results_file,
        open(output_dir / "annotations.tsv", "w", newline="") as annotations_file,
        open(output_dir / "motifs.tsv", "w", newline="") as motifs_file,
        open(output_dir / "domains.tsv", "w", newline="") as domains_file,
        open(output_dir / "nbarc.fa", "w") as nbarc_fasta,
        open(classified_fasta_path, "w") as classified_fasta,
    ):
        results = csv.writer(results_file, delimiter="\t")
        annotations = csv.writer(annotations_file, delimiter="\t")
        motifs = csv.writer(motifs_file, delimiter="\t")
        domains = csv.writer(domains_file, delimiter="\t")

        if command == "nlr":
            results.writerow(
                [
                    "Sequence",
                    "Length",
                    "LRR_Length",
                    "Motifs",
                    "Domains",
                    "Classification",
                    "NBARC_motifs",
                    "MADA",
                    "CJID",
                ]
            )
        elif command == "prr":
            results.writerow(
                [
                    "Sequence",
                    "Length",
                    "Extracellular_Length",
                    "LRR_Length",
                    "Type",
                    "Classification",
                    "Signal_peptide",
                ]
            )

        annotations.writerow(
            ["Sequence", "Domain", "Start", "End", "Accession", "Score", "Source"]
        )
        motifs.writerow(
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
        domains.writerow(["Sequence", "Domain", "Start", "End"])

        for protein in proteins.values():
            if command == "nlr":
                results.writerow(
                    [
                        protein.id,
                        protein.length,
                        protein.lrr_length,
                        protein.motif_string,
                        protein.domain_string,
                        protein.classification,
                        protein.nbarc_motif_count,
                        protein.has_annotation("MADA"),
                        protein.has_annotation("C-JID"),
                    ]
                )
            elif command == "prr":
                results.writerow(
                    [
                        protein.id,
                        protein.length,
                        protein.extracellular_length,
                        protein.lrr_length,
                        protein.type,
                        protein.classification,
                        protein.has_annotation("signal_peptide"),
                    ]
                )

            for annotation in protein.annotations:
                if annotation.type == "domain" and annotation.source != "merged":
                    annotations.writerow(
                        [
                            protein.id,
                            annotation.name,
                            annotation.start,
                            annotation.end,
                            annotation.accession,
                            round(annotation.score, 2)
                            if isinstance(annotation.score, float)
                            else annotation.score,
                            annotation.source,
                        ]
                    )
                elif annotation.type == "domain" and annotation.source == "merged":
                    domains.writerow(
                        [
                            protein.id,
                            annotation.name,
                            annotation.start,
                            annotation.end,
                        ]
                    )
                elif annotation.type == "motif":
                    start, end = annotation.start, annotation.end
                    motifs.writerow(
                        [
                            protein.id,
                            annotation.name,
                            start,
                            round(annotation.score, 2)
                            if isinstance(annotation.score, float)
                            else annotation.score,
                            protein.sequence[max(0, start - 6) : start - 1],
                            protein.sequence[start - 1 : end],
                            protein.sequence[end : end + 5],
                        ]
                    )

            nbarc_domains = [
                a
                for a in protein.annotations
                if a.name == "NBARC" and a.type == "domain" and a.source == "merged"
            ]
            for domain in nbarc_domains:
                nbarc_fasta.write(f">{protein.id}_{domain.start}_{domain.end}\n")
                nbarc_fasta.write(
                    protein.sequence[domain.start - 1 : domain.end] + "\n"
                )

            if protein.classification is not None:
                classified_fasta.write(f">{protein.id}\n{protein.sequence}\n")

            if draw:
                sanitised_id = re.sub(r"[\\/*?:<>|]", "", protein.id)
                draw_svg(
                    protein,
                    output_dir / "plots" / f"{sanitised_id}.svg",
                    command=command,
                )

    if draw:
        plots_dir = output_dir / "plots"
        tar_path = output_dir / "plots.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(plots_dir, arcname="plots")
        shutil.rmtree(plots_dir)
