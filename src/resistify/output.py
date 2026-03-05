import csv
import logging
import re
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
    "RNSB-A": "#FE6100",
    "Walker-B": "#FE6100",
    "RNSB-B": "#FE6100",
    "RNSB-C": "#FE6100",
    "RNSB-D": "#FE6100",
    "GLPL": "#FE6100",
    "MHD": "#FE6100",
    "LxxLxL": "#FFB000",
}


def draw_svg(protein: Protein, output_path: Path):
    elements: list[Element]

    elements = [
        svg.Rect(
            x=0,
            y=30,
            width=protein.length,
            height=40,
            fill="lightgrey",
            stroke="black",
            stroke_width=2,
        ),
    ]

    motifs = [a for a in protein.annotations if a.type == "motif"]
    for motif in motifs:
        elements.append(
            svg.Line(
                x1=motif.start,
                x2=motif.start,
                y1=10,
                y2=30,
                stroke="black",
                stroke_width=2,
            ),
        )
        elements.append(
            svg.Text(
                x=motif.start,
                y=20,
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
    for domain in domains:
        elements.append(
            svg.Rect(
                x=domain.start,
                y=30,
                width=domain.end - domain.start + 1,
                height=40,
                fill=COLOUR_PALETTE.get(domain.name, "grey"),
                stroke="black",
                stroke_width=2,
            )
        )
        elements.append(
            svg.Text(
                x=domain.start + (domain.end - domain.start) / 2,
                y=50,
                text=domain.name,
                font_size=24,
                font_family="sans-serif",
                text_anchor="middle",
                dominant_baseline="middle",
            )
        )

    canvas = svg.SVG(width=protein.length, height=100, elements=elements)

    with open(output_path, "w") as f:
        f.write(str(canvas))


def save_results(proteins: dict[str, Protein], output_dir: Path, command: str):
    logger.info("Saving results")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "plots").mkdir(parents=True, exist_ok=True)

    with (
        open(output_dir / "results.tsv", "w", newline="") as results_file,
        open(output_dir / "annotations.tsv", "w", newline="") as annotations_file,
        open(output_dir / "motifs.tsv", "w", newline="") as motifs_file,
        open(output_dir / "domains.tsv", "w", newline="") as domains_file,
        open(output_dir / "nbarc.fa", "w") as nbarc_fasta,
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
            ["Sequence", "Domain", "Start", "End", "E_value", "Score", "Source"]
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
                            annotation.score,
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
                            annotation.score,
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

            draw_svg(
                protein,
                output_dir / "plots" / f"{re.sub(r'[\\/*?:<>|]', '', protein.id)}.svg",
            )
