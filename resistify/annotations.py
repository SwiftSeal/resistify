import logging

DUPLICATE_GAP = 100

short_IDs = {
    "CC": "C",
    "RPW8": "R",
    "TIR": "T",
    "NB-ARC": "N",
    "LRR": "L",
    "MADA": "m",
    "C-JID": "j",
}


class Sequence:
    def __init__(self, sequence):
        self.sequence = sequence
        self.classification = None
        self.mada = False
        self.cjid = False
        self.domain_string = ""
        self.annotations = []
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
            "bDaD1": []
        }

    def add_annotation(self, annotation):
        self.annotations.append(annotation)

    def add_motif(self, motif):
        self.motifs[motif.classification].append(motif)

    def merge_annotations(self):
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
                if annotation.start <= merged_sublist[-1].end + DUPLICATE_GAP:
                    merged_sublist[-1].end = annotation.end
                else:
                    merged_sublist.append(annotation)
            merged_annotations += merged_sublist

        # sort merged annotations by start position
        merged_annotations.sort(key=lambda x: x.start)
        self.annotations = merged_annotations

    def classify(self):
        domain_string = ""
        for annotation in self.annotations:
            domain_string += short_IDs[annotation.domain]

        # classify based on primary architecture
        if "CN" in domain_string:
            self.classification = "CN"
        elif "RN" in domain_string:
            self.classification = "RN"
        elif "TN" in domain_string:
            self.classification = "TN"
        elif "N" in domain_string:
            self.classification = "N"
        else:
            self.classification = None

    def reclassify(self, lrr_gap, lrr_length):
        """
        Reclassify with new LRR annotations.
        """
        # Add CC annotation from motif if no N terminal annotation
        if self.classification == "N":
            for annotation in self.annotations:
                if annotation.domain == "NB-ARC":
                    nbarc_start = annotation.start
                    break
            for motif in self.motifs["extEDVID"]:
                if motif.position < nbarc_start:
                    self.add_annotation(
                        Annotation("CC", motif.position, motif.position + 1)
                    )
            TIR_motif_IDs = ["aA", "aC", "aD3", "bA", "bC", "bDaD1"]
            TIR_motifs = [item for motif in TIR_motif_IDs for item in self.motifs[motif] if item.position < nbarc_start]
            # TIR motifs are pretty conserved, seems okay to take 1 as sufficient evidence
            if len(TIR_motifs) > 0:
                TIR_motifs.sort(key=lambda x: x.position)
                self.add_annotation(
                    Annotation("TIR", TIR_motifs[0].position, TIR_motifs[-1].position)
                )


        # Add LRR annotation if there are more than 3 LRR motifs
        if len(self.motifs["LxxLxL"]) > lrr_length - 1:
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
                        self.add_annotation(Annotation("LRR", start, end))
                    start = motif.position
                    end = motif.position
                    count = 0

            if count >= 3:
                self.add_annotation(Annotation("LRR", start, end))

        sorted_annotations = sorted(self.annotations, key=lambda x: x.start)

        domain_string = ""
        for annotation in sorted_annotations:
            domain_string += short_IDs[annotation.domain]

        self.domain_string = domain_string

        # check for MADA and C-JID
        if "m" in domain_string:
            self.mada = True
        if "j" in domain_string:
            self.cjid = True

        # classify based on primary architecture
        # Does order matter?
        if "CNL" in domain_string:
            self.classification = "CNL"
        elif "RNL" in domain_string:
            self.classification = "RNL"
        elif "TNL" in domain_string:
            self.classification = "TNL"
        elif "NL" in domain_string:
            self.classification = "NL"
        elif "CN" in domain_string:
            self.classification = "CN"
        elif "TN" in domain_string:
            self.classification = "TN"
        else:
            return


class Annotation:
    def __init__(self, domain, start, end):
        self.domain = domain
        self.start = start
        self.end = end


class Motif:
    classification = str
    probability = float
    position = int

    def __init__(self, classification, probability, position):
        self.classification = classification
        self.probability = probability
        self.position = int(position)
