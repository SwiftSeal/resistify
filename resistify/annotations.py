from __future__ import annotations
from resistify._loguru import logger
from dataclasses import dataclass, field

nlr_classifications = ["RNL", "CNL", "TNL", "RN", "CN", "TN", "NL", "N"]

rlp_external_domains = [
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

short_IDs = {
    "CC": "C",
    "RPW8": "R",
    "TIR": "T",
    "NB-ARC": "N",
    "LRR": "L",
    "MADA": "m",
    "C-JID": "j",
}

# Domains that should be merged if adjacent
# These will also be reported in the domain table
DOMAINS_TO_MERGE = [
    "CC",
    "RPW8",
    "TIR",
    "NB-ARC",
    "LRR",
    "MADA",
    "C-JID",
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
    "PKinase",
    "alpha_inwards",  # don't care about the other topology domains
    "signal_peptide",
]

MOTIF_TRANSLATION = {
    "extEDVID": "C",
    "bA": "T",
    "aA": "T",
    "bC": "T",
    "aC": "T",
    "bDaD1": "T",
    "aD3": "T",
    "VG": "N",
    "P-loop": "N",
    "RNSB-A": "N",
    "Walker-B": "N",
    "RNSB-B": "N",
    "RNSB-C": "N",
    "RNSB-D": "N",
    "GLPL": "N",
    "MHD": "N",
    "LxxLxL": "L",
}


@dataclass
class Annotation:
    domain: str
    start: int
    end: int
    evalue: float | None = None
    score: float | None = None
    source: str | None = None


@dataclass
class Motif:
    type: str
    probability: float
    position: int
    parent_sequence: Sequence
    length: int

    @property
    def domain(self):
        """
        Return the domain type of the motif.
        """
        return MOTIF_TRANSLATION.get(self.type, "?")

    @property
    def sequence(self):
        return self.parent_sequence.seq[self.position : self.position + self.length]

    @property
    def upstream(self, flank=5):
        return self.parent_sequence.seq[max(0, self.position - flank) : self.position]

    @property
    def downstream(self, flank=5):
        end = self.position + self.length
        return self.parent_sequence.seq[end : end + flank]


@dataclass
class Sequence:
    id: str
    seq: str
    type: str | None = None
    classification: str | None = None
    annotations: list[Annotation] = field(default_factory=list)
    merged_annotations: list[Annotation] = field(default_factory=list)
    motifs: list[Motif] = field(default_factory=list)
    cc_probs: list[float] = field(default_factory=list)
    transmembrane_predictions: str | None = None

    @property
    def motif_string(self):
        motif_string = ""
        for motif in self.motifs:
            motif_string += motif.domain
        return motif_string

    @property
    def domain_string(self):
        domain_string = ""
        for annotation in self.merged_annotations:
            domain_string += short_IDs[annotation.domain]
        return domain_string

    @property
    def lrr_length(self):
        """
        Calculate the total length of the LRR domains in the sequence.
        Do this for the merged annotations to acount for gap merging.
        """
        length = 0
        for annotation in self.merged_annotations:
            if annotation.domain == "LRR":
                length += annotation.end - annotation.start + 1
        return length

    @property
    def nterminal_sequence(self):
        """
        Return the N-terminal sequence up to the start of the NB-ARC domain.
        Only the first NB-ARC domain is considered.
        """
        for annotation in self.annotations:
            if annotation.domain == "NB-ARC":
                return self.seq[: annotation.start]
        # If no NB-ARC domain, return the whole sequence
        return None

    @property
    def has_nbarc(self):
        if self.annotations is None:
            return False
        else:
            for annotation in self.annotations:
                if annotation.domain == "NB-ARC":
                    return True

    @property
    def has_mada(self):
        for annotation in self.annotations:
            if annotation.domain == "MADA" and annotation.score >= 20:
                return True
        return False

    @property
    def has_madal(self):
        for annotation in self.annotations:
            if annotation.domain == "MADA" and annotation.score < 20:
                return True
        return False

    @property
    def has_cjid(self):
        for annotation in self.annotations:
            if annotation.domain == "C-JID":
                return True
        return False

    @property
    def has_signal_peptide(self):
        for annotation in self.annotations:
            if annotation.domain == "signal_peptide":
                return True
        return False

    @property
    def extracellular_length(self):
        """
        Calculate the length of the extracellular domain of an RLP/RLK.
        """
        if self.type == "RLP" or self.type == "RLK":
            for annotation in self.annotations:
                if annotation.domain == "alpha_inwards":
                    return annotation.start
        else:
            return None

    def add_annotation(
        self,
        domain: str,
        source: str,
        start: int,
        end: int,
        evalue: float | None = None,
        score: float | None = None,
    ):
        """
        Add annotation and keep sorted by start position.
        """
        if start > end:
            logger.error(f"Invalid annotation coordinates for {self.id}")
            return
        logger.debug(f"Adding annotation {domain} to {self.id} from {start} to {end}")
        self.annotations.append(Annotation(domain, start, end, evalue, score, source))
        self.annotations.sort(key=lambda x: x.start)

    def add_motif(self, type, probability, position, length):
        """
        Add a motif to the sequence.
        Motifs are guaranteed to be sorted by position.
        """
        logger.debug(
            f"Adding motif {type} to {self.id} at position {position} with probability {probability}"
        )
        self.motifs.append(Motif(type, probability, position, self, length))
        self.motifs.sort(key=lambda x: x.position)

    def get_motifs(self, domain=None):
        if domain:
            return [m for m in self.motifs if m.domain == domain]
        return self.motifs

    def identify_cc_domains(self):
        """
        Identify CC domains based on Coconat CC probabilites.
        Use a sliding window of size N and predict as CC if mean probability less than X
        """
        window_size = 5
        threshold = 0.7
        start = None

        for i in range(len(self.cc_probs) - window_size + 1):
            # Get the current window
            current_window = self.cc_probs[i : i + window_size]

            # Calculate the average probability for the current window
            avg_prob = sum(current_window) / window_size

            # Check if the average is below the threshold
            if avg_prob > threshold:
                if start is None:
                    start = i  # Mark the start of the dipping region
                end = i + window_size - 1  # Extend the end of the region
            else:
                # If we were in a dipping region and now the condition is false, record the region
                if start is not None:
                    self.add_annotation("CC", "coconat", start, end)
                    start = None  # Reset start for the next region

        # If we ended in a dip region, capture the final one
        if start is not None:
            self.add_annotation("CC", "coconat", start, end)

    def identify_lrr_domains(self, lrr_gap: int, lrr_length: int):
        """
        Identify LRR domains based on LxxLxL motifs.
        """

        lrr_motifs = self.get_motifs("L")

        if len(lrr_motifs) == 0:
            logger.debug(f"No LRR motifs found in {self.id}")
            return

        current_motif = lrr_motifs[0]
        start = current_motif.position
        end = current_motif.position
        count = 0
        for motif in lrr_motifs[1:]:
            if motif.position - end < lrr_gap:
                end = motif.position
                count += 1
            else:
                if count >= lrr_length:
                    self.add_annotation("LRR", "nlrexpress", start, end)
                start = motif.position
                end = motif.position
                count = 0

        if count >= lrr_length:
            self.add_annotation("LRR", "nlrexpress", start, end)

    def classify_nlr(self):
        # create a simplified domain string
        domain_string = self.domain_string

        # collapse adjacent identical domains for classification
        collapsed_domain_string = ""
        if len(domain_string) > 0:
            collapsed_domain_string = [domain_string[0]]
            for domain in domain_string[1:]:
                if domain != collapsed_domain_string[-1]:
                    collapsed_domain_string.append(domain)
            collapsed_domain_string = "".join(collapsed_domain_string)

        logger.debug(
            f"Collapsed domain string for {self.id}: {collapsed_domain_string}"
        )

        # Absolutely mawkit, but catch RC collapsed string which will occur when coconat is applied to rpw8
        collapsed_domain_string = collapsed_domain_string.replace("RC", "R")

        # classify based on primary architecture - first match wins (go team CNL!)
        for classification in nlr_classifications:
            if classification in collapsed_domain_string:
                self.classification = classification
                self.type = "NLR"
                break

        # scavenge for missed classifications with motifs
        # this is all very janky
        if self.classification == "N" or self.classification == "NL":
            # get the start of the NB-ARC domain
            for annotation in self.annotations:
                if annotation.domain == "NB-ARC":
                    nbarc_start = annotation.start
                    break

            CC_motifs = [m for m in self.get_motifs("C") if m.position < nbarc_start]
            TIR_motifs = [m for m in self.get_motifs("T") if m.position < nbarc_start]

            if len(CC_motifs) > 0:
                self.add_annotation(
                    "CC",
                    "nlrexpress",
                    CC_motifs[0].position,
                    CC_motifs[-1].position,
                )
                self.classification = "C" + self.classification
            elif len(TIR_motifs) > 0:
                TIR_motifs.sort(key=lambda x: x.position)
                self.add_annotation(
                    "TIR",
                    "nlrexpress",
                    TIR_motifs[0].position,
                    TIR_motifs[-1].position,
                )
                self.classification = "T" + self.classification

    def is_rlp(self, extracellular_length=50):
        tm_detected = False
        n_terminal_length = 0

        for annotation in self.annotations:
            # Immediately skip beta barrels or outward alpha helices
            if annotation.domain in ["beta_inwards", "beta_outwards", "alpha_outwards"]:
                return False
            elif annotation.domain == "alpha_inwards":
                # If already detected an alpha helix, not single-pass
                if tm_detected:
                    return False
                tm_detected = True
                n_terminal_length = annotation.start

        if not tm_detected:
            return False

        # Use extracellular length threshold to filter out non-RLPs
        if n_terminal_length < extracellular_length:
            return False
        else:
            self.type = "RLP"
            return True

    def classify_rlp(self):
        """
        Extract features indicative of a RLK/RLP based on TMBed topology.
        Update protein with relevant topology information
        """
        for annotation in self.annotations:
            if annotation.domain == "alpha_inwards":
                tm_end = annotation.end
                tm_start = annotation.start

        external_domains = set()
        for annotation in self.annotations:
            if annotation.domain == "PKinase" and annotation.start > tm_end:
                self.type = "RLK"
            if annotation.domain in rlp_external_domains and annotation.end < tm_start:
                external_domains.add(annotation.domain)

        external_domains = ";".join(external_domains)

        # If no external domains identified, set classification to None
        if not external_domains:
            self.classification = None
        else:
            self.classification = external_domains

    def merge_annotations(self, duplicate_gap):
        """
        Merge overlapping annotations of the same domain.
        """
        merged_annotations = []

        for domain in DOMAINS_TO_MERGE:
            # get all annotations of this domain
            domain_sublist = [
                annotation
                for annotation in self.annotations
                if annotation.domain == domain
            ]
            # skip if there are no annotations of this domain
            if len(domain_sublist) == 0:
                continue
            domain_sublist.sort(key=lambda x: x.start)
            # merge overlapping annotations
            merged_sublist = [domain_sublist[0]]
            for annotation in domain_sublist[1:]:
                if annotation.start <= merged_sublist[-1].end + duplicate_gap:
                    merged_sublist[-1].end = annotation.end
                else:
                    merged_sublist.append(annotation)
            merged_annotations += merged_sublist

        # sort merged annotations by start position
        merged_annotations.sort(key=lambda x: x.start)
        self.merged_annotations = merged_annotations
