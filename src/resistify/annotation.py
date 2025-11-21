from dataclasses import dataclass, field
import csv
import logging
from pathlib import Path
import svg
from svg import Element

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
}


@dataclass
class Annotation:
    name: str  # Primary name of domain
    start: int  # 1-based
    end: int  # funnily enough also 1-based
    type: str = "domain"  # can also be motif
    source: str | None = None
    accession: str | None = None
    score: float | None = None

    def __post_init__(self):
        if not (1 <= self.start <= self.end):
            raise ValueError("Domain start must be <= end and both must be positive.")


@dataclass
class Protein:
    id: str
    sequence: str
    classification: str | None = None
    annotations: list[Annotation] = field(default_factory=list)
    cc_probs: list[float] | None = None

    def __post_init__(self):
        self.length = len(self.sequence)

    def get_annotation_by_name(self, annotation_name: str) -> list[Annotation]:
        return [
            annotation
            for annotation in self.annotations
            if annotation.name == annotation_name
        ]

    def add_annotation(self, annotation: Annotation):
        if not (1 <= annotation.start <= annotation.end <= self.length):
            raise ValueError(
                f"Annotation boundaries ({annotation.start}-{annotation.end}) out of bounds for sequence of length {self.length}."
            )
        self.annotations.append(annotation)

    def has_annotation(self, annotation_name: str) -> bool:
        for annotation in self.annotations:
            if annotation.name == annotation_name:
                return True
        return False

    def annotate_lrr(self, lrr_gap: int = 75, lrr_length: int = 4):
        """
        Annotate LRR domains based on existing LxxLxL motifs.
        """
        motif_positions = [a.start - 1 for a in self.annotations if a.name == "LxxLxL"]
        if motif_positions:
            start, end, count = 0, 0, 0
            for pos in motif_positions:
                if pos - end < lrr_gap:
                    # Extend if still within gap limit
                    end = pos
                    count += 1
                else:
                    if count >= lrr_length:
                        # Annotate if big enough
                        self.add_annotation(Annotation("LRR", start + 1, end))
                    start = pos
                    end = pos
                    count = 1
            if count >= lrr_length:
                # Annotate if big enough
                self.add_annotation(Annotation("LRR", start + 1, end))

    def annotate_cc(self):
        """
        Identify CC domains based on Coconat CC probabilites.
        Use a sliding window of size N and predict as CC if mean probability less than X
        """
        window_size: int = 5
        threshold: float = 0.3
        start: None | int = None
        end: None | int = None

        if self.cc_probs is None:
            return

        for i in range(len(self.cc_probs) - window_size + 1):
            # Get the current window
            current_window = self.cc_probs[i : i + window_size]

            # Calculate the average probability for the current window
            avg_prob = sum(current_window) / window_size

            # Check if the average is above the threshold
            if avg_prob > threshold:
                if start is None:
                    start = i  # Mark the start of the region
                end = i + window_size - 1  # Extend the end of the region
            else:
                # If we were in a region and now the condition is false, record the region
                if start is not None:
                    self.add_annotation(
                        Annotation(
                            "CC", start + 1, end + 1, type="domain", source="coconat"
                        )
                    )
                    start = None  # Reset start for the next region

        # If we ended in a region, capture the final one
        if start is not None:
            self.add_annotation(
                Annotation("CC", start + 1, end + 1, type="domain", source="coconat")
            )

    def merge_domains(self, gap: int = 50):
        """
        Merge overlapping or adjacent domain annotations of the same type.
        Add as new annotations.
        """
        domain_types = set(a.name for a in self.annotations if a.type == "domain")
        for domain_type in domain_types:
            domain_annotations = sorted(
                (a for a in self.annotations if a.name == domain_type),
                key=lambda x: x.start,
            )
            current_start, current_end = (
                domain_annotations[0].start,
                domain_annotations[0].end,
            )
            for annotation in domain_annotations[1:]:
                if annotation.start <= current_end + gap:
                    # Overlapping or close enough to merge
                    current_end = max(current_end, annotation.end)
                else:
                    # No overlap, save current and move to next
                    self.add_annotation(
                        Annotation(
                            name=domain_type,
                            start=current_start,
                            end=current_end,
                            type="domain",
                            source="merged",
                        )
                    )
                    current_start, current_end = annotation.start, annotation.end
            # Add the last merged domain
            self.add_annotation(
                Annotation(
                    name=domain_type,
                    start=current_start,
                    end=current_end,
                    type="domain",
                    source="merged",
                )
            )

    def classify_nlr(self):
        """
        Uses available annotation data to classify NLRs according to structure.
        """
        present_domains = {a.name for a in self.annotations if a.type == "domain"}
        if "NBARC" in present_domains:
            if "LRR" in present_domains:
                if "RPW8" in present_domains:
                    self.classification = "RNL"
                elif "CC" in present_domains:
                    self.classification = "CNL"
                elif "TIR" in present_domains:
                    self.classification = "TNL"
                else:
                    self.classification = "NL"
            else:
                if "RPW8" in present_domains:
                    self.classification = "RN"
                elif "CC" in present_domains:
                    self.classification = "CN"
                elif "TIR" in present_domains:
                    self.classification = "TN"
                else:
                    self.classification = "N"
        else:
            # No NB-ARC domain, classify as non-NLR
            self.classification = "non-NLR"

    def draw_svg(self, output_path: Path):
        elements: list[Element]

        elements = [
            svg.Line(
                x1=0, x2=self.length, y1=50, y2=50, stroke="black", stroke_width=2
            ),
        ]

        motifs = [a for a in self.annotations if a.type == "motif"]
        for motif in motifs:
            elements.append(
                svg.Line(
                    x1=motif.start,
                    x2=motif.start,
                    y1=50,
                    y2=20,
                    stroke="black",
                    stroke_width=1,
                ),
            )
            elements.append(
                svg.Text(
                    x=motif.start,
                    y=20,
                    text=motif.name,
                    transform=[svg.Rotate(45, motif.start, 20)],
                    font_size=8,
                    font_family="sans-serif",
                )
            )

        domains = [
            a for a in self.annotations if a.type == "domain" and a.source == "merged"
        ]
        for domain in domains:
            elements.append(
                svg.Rect(
                    x=domain.start,
                    y=30,
                    width=domain.end - domain.start + 1,
                    height=40,
                    fill=COLOUR_PALETTE.get(domain.name, "grey"),
                )
            )
            elements.append(
                svg.Text(
                    x=domain.start + (domain.end - domain.start) / 2,
                    y=50,
                    text=domain.name,
                    font_size=12,
                    font_family="sans-serif",
                    font_weight="bold",
                    text_anchor="middle",
                )
            )

        canvas = svg.SVG(width=self.length, height=100, elements=elements)

        with open(output_path, "w") as f:
            f.write(str(canvas))


def save_results(proteins: dict[str, Protein], output_dir: Path):
    logger.info("Saving results")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "plots").mkdir(parents=True, exist_ok=True)
    results = csv.writer(open(output_dir / "results.tsv", "w"), delimiter="\t")
    annotations = csv.writer(open(output_dir / "annotations.tsv", "w"), delimiter="\t")

    results.writerow(["id", "classification", "length", "num_lrr_motifs"])
    annotations.writerow(["id", "name", "start", "end", "type", "accession", "score"])

    for protein in proteins.values():
        results.writerow(
            [
                protein.id,
                protein.classification,
                protein.length,
                sum(1 for a in protein.annotations if a.name == "LxxLxL"),
            ]
        )
        for annotation in protein.annotations:
            annotations.writerow(
                [
                    protein.id,
                    annotation.name,
                    annotation.start,
                    annotation.end,
                    annotation.type,
                    annotation.source,
                    annotation.accession,
                    annotation.score,
                ]
            )
        protein.draw_svg(output_dir / "plots" / f"{protein.id}.svg")
