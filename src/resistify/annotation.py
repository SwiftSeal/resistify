from dataclasses import dataclass, field
import csv
from pathlib import Path
from resistify.console import console


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

    @property
    def has_annotation(self, annotation_name: str) -> bool:
        for annotation in self.annotations:
            if annotation.name == annotation_name:
                return True
        return False


def classify_nlrs(proteins: list[Protein]):
    """
    Uses available annotation data to classify NLRs according to structure.
    """
    console.log("Classifying NLRs")
    for protein in proteins:
        present_domains = {a.name for a in protein.annotations}
        if "NBARC" in present_domains:
            if "LRR" in present_domains:
                if "RPW8" in present_domains:
                    protein.classification = "RNL"
                elif "CC" in present_domains:
                    protein.classification = "CNL"
                elif "TIR" in present_domains:
                    protein.classification = "TNL"
                else:
                    protein.classification = "NL"
            else:
                if "RPW8" in present_domains:
                    protein.classification = "RN"
                elif "CC" in present_domains:
                    protein.classification = "CN"
                elif "TIR" in present_domains:
                    protein.classification = "TN"
                else:
                    protein.classification = "N"
        else:
            # No NB-ARC domain, classify as non-NLR
            protein.classification = "non-NLR"


def save_results(proteins: list[Protein], output_dir: Path):
    console.log("Saving results")
    results = csv.writer(open(output_dir / "results.tsv", "w"), delimiter="\t")
    annotations = csv.writer(open(output_dir / "annotations.tsv", "w"), delimiter="\t")

    results.writerow(["id", "classification", "length", "num_lrr_motifs"])
    annotations.writerow(["id", "name", "start", "end", "type", "accession", "score"])

    for protein in proteins:
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
