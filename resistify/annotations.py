import logging

log = logging.getLogger(__name__)

short_IDs = {
    "CC": "C",
    "RPW8": "R",
    "TIR": "T",
    "NB-ARC": "N",
    "LRR": "L",
    "MADA": "m",
    "C-JID": "j",
}

motif_translation = {
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


class Sequence:
    def __init__(self, id, sequence):
        self.id = id
        self.sequence = sequence
        self.classification = None
        self.mada = False
        self.madal = False
        self.cjid = False
        self.domain_string = ""
        self.annotations = []
        self.merged_annotations = []
        self.motifs = {
            "extEDVID": [],
            "VG": [],
            "P-loop": [],
            "RNSB-A": [],
            "Walker-B": [],
            "RNSB-B": [],
            "RNSB-C": [],
            "RNSB-D": [],
            "GLPL": [],
            "MHD": [],
            "LxxLxL": [],
            "aA": [],
            "aC": [],
            "aD3": [],
            "bA": [],
            "bC": [],
            "bDaD1": [],
        }
        self.cc_probs = []

    def motif_string(self):
        sorted_motifs = [item for sublist in self.motifs.values() for item in sublist]
        sorted_motifs.sort(key=lambda x: x.position)
        motif_string = ""
        for motif in sorted_motifs:
            motif_string += motif_translation[motif.classification]
        return motif_string

    def add_annotation(self, annotation):
        """
        Add annotation and keep sorted by start position.
        """
        self.annotations.append(annotation)
        self.annotations.sort(key=lambda x: x.start)

    def add_motif(self, motif):
        self.motifs[motif.classification].append(motif)

    def identify_cc_domains(self):
        """
        Identify CC domains based on Coconat CC probabilites.
        Use a sliding window of size N and predict as CC if mean probability less than X
        """
        window_size = 3
        threshold = 0.7
        start = None

        for i in range(len(self.cc_probs) - window_size + 1):
            # Get the current window
            current_window = self.cc_probs[i : i + window_size]

            # Calculate the average probability for the current window
            avg_prob = sum(current_window) / window_size

            # Check if the average is below the threshold
            if avg_prob < threshold:
                if start is None:
                    start = i  # Mark the start of the dipping region
                end = i + window_size - 1  # Extend the end of the region
            else:
                # If we were in a dipping region and now the condition is false, record the region
                if start is not None:
                    log.debug(f"Adding CC domain in {self.id} from {start} to {end}")
                    self.add_annotation(
                        Annotation("CC", start, end, None, None, "Coconat")
                    )
                    start = None  # Reset start for the next region

        # If we ended in a dip region, capture the final one
        if start is not None:
            log.debug(f"Adding CC domain in {self.id} from {start} to {end}")
            self.add_annotation(Annotation("CC", start, end, None, None, "Coconat"))

    def identify_lrr_domains(self, lrr_gap, lrr_length):
        """
        Identify LRR domains based on LxxLxL motifs.
        """
        if len(self.motifs["LxxLxL"]) == 0:
            return

        sorted_lrr = sorted(self.motifs["LxxLxL"], key=lambda x: x.position)

        current_motif = sorted_lrr[0]
        start = current_motif.position
        end = current_motif.position
        count = 0
        for motif in sorted_lrr[1:]:
            if motif.position - end < lrr_gap:
                end = motif.position
                count += 1
            else:
                if count >= lrr_length:
                    self.add_annotation(
                        Annotation("LRR", start, end, "NA", "NA", "NLRexpress")
                    )
                start = motif.position
                end = motif.position
                count = 0

        if count >= lrr_length:
            self.add_annotation(Annotation("LRR", start, end, "NA", "NA", "NLRexpress"))

    def get_nterminal(self):
        for annotation in self.annotations:
            if annotation.domain == "NB-ARC":
                return self.sequence[: annotation.start]

    def classify(self):
        # create a simplified domain string
        domain_string = ""
        for annotation in self.annotations:
            # skip non-core and flag
            if annotation.domain == "MADA":
                if annotation.score >= 20:
                    self.mada = True
                else:
                    self.madal = True
            elif annotation.domain == "C-JID":
                self.cjid = True
            else:
                domain_string += short_IDs[annotation.domain]

        # collapse adjacent identical domains for classification
        collapsed_domain_string = ""
        if len(domain_string) > 0:
            collapsed_domain_string = [domain_string[0]]
            for domain in domain_string[1:]:
                if domain != collapsed_domain_string[-1]:
                    collapsed_domain_string.append(domain)
            collapsed_domain_string = "".join(collapsed_domain_string)

        log.debug(f"Collapsed domain string for {self.id}: {collapsed_domain_string}")

        # classify based on primary architecture - first match wins (go team CNL!)
        classifications = ["RNL", "CNL", "TNL", "RN", "CN", "TN", "NL", "N"]
        for classification in classifications:
            if classification in collapsed_domain_string:
                self.classification = classification
                break

        # scavenge for missed classifications with motifs
        # this is all very janky
        if self.classification == "N" or self.classification == "NL":
            # get the start of the NB-ARC domain
            for annotation in self.annotations:
                if annotation.domain == "NB-ARC":
                    nbarc_start = annotation.start
                    break
            for motif in self.motifs["extEDVID"]:
                if motif.position < nbarc_start:
                    self.add_annotation(
                        Annotation(
                            "CC",
                            motif.position,
                            motif.position + 1,
                            None,
                            None,
                            "NLRexpress",
                        )
                    )
                    self.classification = "C" + self.classification
                    continue
            TIR_motif_IDs = ["aA", "aC", "aD3", "bA", "bC", "bDaD1"]
            TIR_motifs = [
                item
                for motif in TIR_motif_IDs
                for item in self.motifs[motif]
                if item.position < nbarc_start
            ]
            # TIR motifs are pretty conserved, seems okay to take 1 as sufficient evidence
            if len(TIR_motifs) > 0:
                TIR_motifs.sort(key=lambda x: x.position)
                self.add_annotation(
                    Annotation(
                        "TIR",
                        TIR_motifs[0].position,
                        TIR_motifs[-1].position,
                        None,
                        None,
                        "NLRexpress",
                    )
                )
                self.classification = "T" + self.classification

    def merge_annotations(self, duplicate_gap):
        """
        Merge overlapping annotations of the same domain.
        Don't trust e-values etc - they're inherited from the first annotation.
        """
        merged_annotations = []
        for domain in short_IDs:
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

        domain_string = ""
        for annotation in merged_annotations:
            domain_string += short_IDs[annotation.domain]
        self.domain_string = domain_string


class Annotation:
    def __init__(self, domain, start, end, evalue, score, source):
        self.domain = domain
        self.start = start
        self.end = end
        self.evalue = evalue
        self.score = score
        self.source = source


class Motif:
    classification = str
    probability = float
    position = int

    def __init__(self, classification, probability, position):
        self.classification = classification
        self.probability = probability
        self.position = int(position)
