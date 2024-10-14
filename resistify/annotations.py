import logging

log = logging.getLogger(__name__)

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
    def __init__(self, sequence):
        self.sequence = sequence
        self.classification = None
        self.mada = False
        self.madal = False
        self.cjid = False
        self.domain_string = ""
        self.motif_string = ""
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
            "bDaD1": [],
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
            # skip non-core
            if annotation.domain == "MADA":
                if annotation.score >= 20:
                    self.mada = True
                else:
                    self.madal = True
            elif annotation.domain == "C-JID":
                self.cjid = True
            else:
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
        # Create sorted motif string
        sorted_motifs = [item for sublist in self.motifs.values() for item in sublist]
        sorted_motifs.sort(key=lambda x: x.position)

        # write motif string
        for motif in sorted_motifs:
            self.motif_string += motif_translation[motif.classification]

        # Add CC annotation from motif if no N terminal annotation
        if self.classification == "N":
            for annotation in self.annotations:
                if annotation.domain == "NB-ARC":
                    nbarc_start = annotation.start
                    break
            for motif in self.motifs["extEDVID"]:
                if motif.position < nbarc_start:
                    self.add_annotation(
                        Annotation("CC", motif.position, motif.position + 1, "NA", "NA")
                    )
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
                        "NA",
                        "NA",
                    )
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
                        self.add_annotation(Annotation("LRR", start, end, "NA", "NA"))
                    start = motif.position
                    end = motif.position
                    count = 0

            if count >= 3:
                self.add_annotation(Annotation("LRR", start, end, "NA", "NA"))

        sorted_annotations = sorted(self.annotations, key=lambda x: x.start)

        domain_string = ""
        for annotation in sorted_annotations:
            # skip non-core
            if annotation.domain != "MADA" and annotation.domain != "C-JID":
                domain_string += short_IDs[annotation.domain]

        self.domain_string = domain_string

        # collapse adjacent domains
        if len(domain_string) > 0:
            collapsed_domain_string = [domain_string[0]]
            for domain in domain_string[1:]:
                if domain != collapsed_domain_string[-1]:
                    collapsed_domain_string.append(domain)
            domain_string = "".join(collapsed_domain_string)


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
    def __init__(self, domain, start, end, evalue, score):
        self.domain = domain
        self.start = start
        self.end = end
        self.evalue = evalue
        self.score = score


class Motif:
    classification = str
    probability = float
    position = int

    def __init__(self, classification, probability, position):
        self.classification = classification
        self.probability = probability
        self.position = int(position)
