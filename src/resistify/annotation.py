from dataclasses import dataclass, field
import bisect
import logging

logger = logging.getLogger(__name__)


SHORT_MOTIF_IDS = {
    "extEDVID": "C",
    "bA": "T",
    "aA": "T",
    "bC": "T",
    "aC": "T",
    "bDaD1": "T",
    "aD3": "T",
    "VG": "N",
    "P-loop": "N",
    "RNBS-A": "N",
    "Walker-B": "N",
    "RNBS-B": "N",
    "RNBS-C": "N",
    "RNBS-D": "N",
    "GLPL": "N",
    "MHD": "N",
    "LxxLxL": "L",
}

SHORT_DOMAIN_IDS = {
    "CC": "C",
    "RPW8": "R",
    "TIR": "T",
    "NBARC": "N",
    "LRR": "L",
    "MADA": "m",
    "C-JID": "j",
}

RLP_EXTERNAL_DOMAINS = [
    "LRR",
    "G-LecRLK",
    "L-LecRLK",
    "C-LecRLK",
    "WAK",
    "CrRLK1L",
    "LysM",
    "CRK",
    "Thaumatin",
    "CR-like",
    "SPARK",
    "GH18",
    "GH19",
    "CAP",
    "PRIMA1",
]


@dataclass
class Annotation:
    name: str  # Primary name of domain
    start: int  # 1-based
    end: int  # funnily enough also 1-based
    type: str | None = None  # e.g., 'domain' or 'motif'
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
    type: str | None = None  # Can be NLR, RLK, RLP, or None
    classification: str | None = None
    annotations: list[Annotation] = field(default_factory=list)
    cc_probs: list[float] | None = None
    transmembrane_predictions: str | None = None

    def __post_init__(self):
        self.length = len(self.sequence)

    @property
    def motifs(self) -> list[Annotation]:
        """
        Return a list of motif annotations.
        """
        return [a for a in self.annotations if a.type == "motif"]

    @property
    def lrr_length(self) -> int:
        """
        Calculate and return the total length of all LRR domain annotations.
        """
        total_length = 0
        for annotation in self.annotations:
            if annotation.name == "LRR" and annotation.type == "domain":
                total_length += annotation.end - annotation.start + 1
        return total_length

    @property
    def extracellular_length(self) -> int:
        """
        Calculate and return the length of the extracellular region based on transmembrane predictions.
        """
        if self.transmembrane_predictions is None:
            return 0

        tm_end = None
        for annotation in self.annotations:
            if annotation.name == "alpha_inwards":
                tm_end = annotation.start
                break

        if tm_end is None:
            return 0

        return tm_end - 1

    @property
    def nbarc_start(self) -> int | None:
        """
        Return the start position of the NBARC domain, if present.
        """
        for annotation in self.annotations:
            if annotation.name == "NBARC" and annotation.type == "domain":
                return annotation.start
        return None

    @property
    def motif_string(self) -> str:
        """
        Returns a string of single-letter codes representing NLR motifs.
        """
        string = ""
        for motif in self.motifs:
            string += SHORT_MOTIF_IDS.get(motif.name, "")
        return string

    @property
    def domain_string(self) -> str:
        """
        Returns a string of single-letter codes representing NLR domains.
        Uses the merged domain annotations only.
        """
        string = ""
        domain_annotations = [
            a for a in self.annotations if a.type == "domain" and a.source == "merged"
        ]
        for domain in domain_annotations:
            string += SHORT_DOMAIN_IDS.get(domain.name, "")
        return string

    @property
    def nbarc_motif_count(self) -> int:
        """
        Count the number of NBARC motifs present in the protein.
        """
        count = 0
        for motif in self.motifs:
            if motif.name in [
                "VG",
                "P-loop",
                "RNSB-A",
                "Walker-B",
                "RNSB-B",
                "RNSB-C",
                "RNSB-D",
                "GLPL",
                "MHD",
            ]:
                count += 1
        return count

    def add_annotation(self, annotation: Annotation):
        """
        Add an annotation to the protein, ensuring that it is sorted by start position.
        """
        if not (1 <= annotation.start <= annotation.end <= self.length):
            raise ValueError(
                f"Annotation boundaries ({annotation.start}-{annotation.end}) out of bounds for sequence of length {self.length}"
            )
        bisect.insort(self.annotations, annotation, key=lambda x: x.start)

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
                        self.add_annotation(
                            Annotation(
                                "LRR", start + 1, end, type="domain", source="resistify"
                            )
                        )
                    start = pos
                    end = pos
                    count = 1
            if count >= lrr_length:
                # Annotate if big enough
                self.add_annotation(
                    Annotation("LRR", start + 1, end, type="domain", source="resistify")
                )

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

        # Add these last to avoid modifying the list while iterating
        merged_domains = []

        for domain_type in domain_types:
            domain_annotations = [a for a in self.annotations if a.name == domain_type]

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
                    merged_domains.append(
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
            merged_domains.append(
                Annotation(
                    name=domain_type,
                    start=current_start,
                    end=current_end,
                    type="domain",
                    source="merged",
                )
            )

        for merged_domain in merged_domains:
            self.add_annotation(merged_domain)

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
            self.classification = None

    def is_nlr(self):
        """
        Determines if a protein is an NLR based on presence of NBARC domain.
        """
        return self.has_annotation("NBARC")

    def is_rlp(self, minimum_extracellular_length=50):
        tm_detected = False
        n_terminal_length = 0

        for annotation in self.annotations:
            # Immediately skip beta barrels or outward alpha helices
            if annotation.name in ["beta_inwards", "beta_outwards", "alpha_outwards"]:
                return False
            elif annotation.name == "alpha_inwards":
                # If already detected an alpha helix, not single-pass
                if tm_detected:
                    return False
                tm_detected = True
                n_terminal_length = annotation.start

        if not tm_detected:
            return False

        # Use extracellular length threshold to filter out non-RLPs
        if n_terminal_length < minimum_extracellular_length:
            return False
        else:
            return True

    def classify_rlp(self):
        """
        Extract features indicative of a RLK/RLP based on TMBed topology.
        Update protein with relevant topology information
        """
        # Get the positions of the transmembrane domain
        for annotation in self.annotations:
            if annotation.name == "alpha_inwards":
                tm_end = annotation.end
                tm_start = annotation.start

        # Iterate through extracellular annotations and use to set the classification. Motifs and other annotations will be picked up in this, but as long as they are not in the dictionary it doesn't really matter.
        external_domains = set()
        for annotation in self.annotations:
            if annotation.name == "PKinase" and annotation.start > tm_end:
                self.type = "RLK"
            if annotation.name in RLP_EXTERNAL_DOMAINS and annotation.end < tm_start:
                external_domains.add(annotation.name)

        external_domains = ";".join(external_domains)

        # If no external domains identified, set classification to None
        if not external_domains:
            self.classification = None
        else:
            self.classification = external_domains
