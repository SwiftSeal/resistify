short_IDs = {
    "CC": "C",
    "RPW8": "R",
    "TIR": "T",
    "NB-ARC": "N",
}

class Sequence:
    def __init__(self, sequence):
        self.sequence = sequence
        self.classification = None
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
        }


    def add_annotation(self, annotation):
        self.annotations.append(annotation)

    def add_motif(self, motif):
        self.motifs[motif.type].append(motif)
    
    def merge_annotations(self):
        merged_annotations = []
        for domain in short_IDs:
            # get all annotations of this domain
            domain_sublist = [annotation for annotation in self.annotations if annotation.domain == domain]
            # skip if there are no annotations of this domain
            if len(domain_sublist) == 0:
                continue
            domain_sublist.sort(key=lambda x: x.start)
            # merge overlapping annotations
            merged_sublist = [domain_sublist[0]]
            for annotation in domain_sublist[1:]:
                if annotation.start <= merged_sublist[-1].end:
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
        self.position = position