import logging

log = logging.getLogger(__name__)

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
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq
        self.type = None
        self.classification = None
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
        self.transmembrane_predictions = None
        self.signal_peptide = False

    @property
    def motif_string(self):
        sorted_motifs = [item for sublist in self.motifs.values() for item in sublist]
        sorted_motifs.sort(key=lambda x: x.position)
        motif_string = ""
        for motif in sorted_motifs:
            motif_string += motif_translation[motif.classification]
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
        for annotation in self.annotations:
            if annotation.domain == "NB-ARC":
                return True
        return False
    
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

    def add_annotation(self, domain, start, end, evalue, score, source):
        """
        Add annotation and keep sorted by start position.
        """
        self.annotations.append(Annotation(domain, start, end, evalue, score, source))
        self.annotations.sort(key=lambda x: x.start)

    def add_motif(self, predictor, value, i):
        self.motifs[predictor].append(Motif(predictor, value, i))

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
                    log.debug(f"Adding CC domain in {self.id} from {start} to {end}")
                    self.add_annotation("CC", start, end, "NA", "NA", "Coconat")
                    start = None  # Reset start for the next region

        # If we ended in a dip region, capture the final one
        if start is not None:
            log.debug(f"Adding CC domain in {self.id} from {start} to {end}")
            self.add_annotation("CC", start, end, "NA", "NA", "Coconat")

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
                    self.add_annotation("LRR", start, end, "NA", "NA", "NLRexpress")
                start = motif.position
                end = motif.position
                count = 0

        if count >= lrr_length:
            self.add_annotation("LRR", start, end, "NA", "NA", "NLRexpress")

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

        log.debug(f"Collapsed domain string for {self.id}: {collapsed_domain_string}")

        # Absolutely mawkit, but catch RC collapsed string which will occur when coconat is applied to rpw8
        collapsed_domain_string = collapsed_domain_string.replace("RC", "R")

        # classify based on primary architecture - first match wins (go team CNL!)
        for classification in nlr_classifications:
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
                        "CC",
                        motif.position,
                        motif.position + 1,
                        "NA",
                        "NA",
                        "NLRexpress",
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
                    "TIR",
                    TIR_motifs[0].position,
                    TIR_motifs[-1].position,
                    "NA",
                    "NA",
                    "NLRexpress",
                )
                self.classification = "T" + self.classification

        if self.classification in nlr_classifications:
            self.type = "NLR"

    def is_rlp(self):
        tm_detected = False
        inside_count, outside_count = 0, 0
        n_terminal_length = 0

        # If Beta-helixes or IN -> OUT transitions are detected, assume not relevant
        if any(state in self.transmembrane_predictions for state in ["B", "b", "H"]):
            return False

        for i, state in enumerate(self.transmembrane_predictions):
            # set initial states
            if i == 0:
                previous_state = state
                state_start = 0
                continue

            if not tm_detected:
                n_terminal_length += 1
                if state == "i":
                    inside_count += 1
                elif state == "o":
                    outside_count += 1

            if state != previous_state:
                length = i - state_start
                if previous_state == "S" and length > 5:
                    self.signal_peptide = True
                elif previous_state == "h":
                    if tm_detected:
                        return False
                    self.add_annotation(
                        "transmembrane", state_start, i - 1, "NA", "NA", "tmbed"
                    )
                    tm_detected = True

                previous_state = state
                state_start = i

        if previous_state == "h" and not tm_detected:
            self.add_annotation(
                "transmembrane", state_start, i - 1, "NA", "NA", "tmbed"
            )
            tm_detected = True

        if tm_detected is False:
            return False

        if n_terminal_length < 50:
            return False

        if inside_count > outside_count:
            # super rough
            return False

        # As all passed, assume we have a single-pass alpha helix protein
        self.type = "RLP"
        return True

    def classify_rlp(self):
        """
        Extract features indicative of a RLK/RLP based on TMBed topology.
        Update protein with relevant topology information
        """
        for annotation in self.annotations:
            if annotation.domain == "transmembrane":
                tm_end = annotation.end
                tm_start = annotation.start

        external_domains = set()
        for annotation in self.annotations:
            if annotation.domain == "PKinase" and annotation.start > tm_end:
                self.type = "RLK"
            if annotation.domain in rlp_external_domains and annotation.end < tm_start:
                external_domains.add(annotation.domain)

        external_domains = ";".join(external_domains)

        self.classification = external_domains

    def merge_annotations(self, duplicate_gap):
        """
        Merge overlapping annotations of the same domain.
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
