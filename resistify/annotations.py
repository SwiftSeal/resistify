import logging
import csv

classifications = {
    "SSF46785": "Winged-helix",
    "G3DSA:1.10.10.10": "Winged-helix",
    "SSF53474": "a/b",
    "G3DSA:3.40.50.1820": "a/b",
    "PF06760": "MLKL",
    "PF06985": "HET",
    "SM00114": "CARD",
    "PF00619": "CARD",
    "cd08330": "CARD",
    "cd08323": "CARD",
    "cd08788": "CARD",
    "PF18461": "CARD",
    "G3DSA:1.10.533.20": "CARD",
    "PF18461": "CARD",
    "G3DSA:1.10.533.20": "CARD",
    "PF02758": "PYD",
    "PS50824": "PYD",
    "SM01289": "PYD",
    "cd08320": "PYD",
    "cd08321": "PYD",
    "PF00653": "BIR",
    "SM00238": "BIR",
    "PS01282": "BIR",
    "PS50143": "BIR",
    "cd00022": "BIR",
    "PF14484": "FISNA",
    "SM01288": "FISNA",
    "PF17107": "NN",
    "PF17779": "nacht-winged-helix",
    "PF17776": "nacht-winged-helix",
    "PF05729": "NACHT",
    "PS50837": "NACHT",
    "PF12061": "R1",
    "cd14798": "CC",
    "PF18052": "CC",
    "G3DSA:1.20.5.4130": "CC",
    "PS51153": "RPW8",
    "PF05659": "RPW8",
    "G3DSA:1.20.930.20": "CblN",
    "G3DSA:3.40.50.10140": "TIR",
    "PF01582": "TIR",
    "PS50104": "TIR",
    "PF13676": "TIR",
    "SM00255": "TIR",
    "SSF52200": "TIR",
    "G3DSA:3.40.50.300": "PLOOP",
    "SSF52540": "PLOOP",
    "G3DSA:1.10.8.430": "NBARC",
    "PF00931": "NBARC",
    "PR00364": "disease",
    "G3DSA:3.80.10.10": "LRR",
    "PF08263": "LRR",
    "PF07723": "LRR",
    "PF07725": "LRR",
    "PF12799": "LRR",
    "PF13306": "LRR",
    "PF00560": "LRR",
    "PF13516": "LRR",
    "PF13855": "LRR",
    "SSF52047": "LRR",
    "SSF52058": "LRR",
    "SM00367": "LRR",
    "SM00368": "LRR",
    "SM00369": "LRR",
    "PF18837": "LRR",
    "PF01463": "LRR",
    "SM00082": "LRR",
    "SM00013": "LRR",
    "PF01462": "LRR",
    "PF18831": "LRR",
    "PF18805": "LRR",
    "CJID": "CJID",
    "PF20160": "CJID",
    "Motif 2": "rnbs-d",
    "Motif 17": "cc-motif",
    "Motif 16": "cc-motif",
    "Motif 14": "cc-motif",
    "Motif 1": "other-motif",
    "Motif 3": "other-motif",
    "Motif 5": "other-motif",
    "Motif 6": "other-motif",
    "Motif 10": "other-motif",
    "Motif 12": "other-motif",
    "Motif 8": "linker-MHD",
    "Motif 7": "linker-MHD",
    "PS50297": "ANK",
    "PF12796": "ANK",
    "PF11929": "ANK",
    "SM00248": "ANK",
    "PS50088": "ANK",
    "PF00023": "ANK",
    "PF13606": "ANK",
    "G3DSA:1.25.40.20": "ANK",
    "SSF48403": "ANK",
    "G3DSA:2.130.10.10": "WD40",
    "SSF50978": "WD40",
    "PS50294": "WD40",
    "PF16756": "WD40",
    "PF16529": "WD40",
    "PF12894": "WD40",
    "SM00320": "WD40",
    "PF00400": "WD40",
    "PS50082": "WD40",
    "cd00200": "WD40",
    "PS00678": "WD40",
    "SSF48371": "ARM",
    "G3DSA:1.25.10.10": "ARM",
    "SM00185": "ARM",
    "PF00514": "ARM",
    "PS50176": "ARM",
    "PF04826": "ARM",
    "PF02985": "ARM",
    "PF01602": "ARM",
    "PF13646": "ARM",
    "G3DSA:1.25.40.10": "TPR",
    "SSF48452": "TPR",
    "PF00515": "TPR",
    "PF07719": "TPR",
    "PF07720": "TPR",
    "PF07721": "TPR",
    "PF12688": "TPR",
    "PF13374": "TPR",
    "PF13424": "TPR",
    "PF09976": "TPR",
    "SM00028": "TPR",
    "PS50005": "TPR",
    "PF13176": "TPR",
    "PF13181": "TPR",
    "PF13174": "TPR",
    "PS50293": "TPR",
    "SM00671": "TPR",
    "PF08238": "TPR",
    "PF17874": "TPR",
    "PF18391": "TPR",
    "PF10516": "TPR",
    "PF18768": "TPR",
    "G3DSA:1.20.58.320": "TPR",
    "SM00090": "PKin",
    "cd05144": "PKin",
    "PF00069": "PKin",
    "PS50011": "PKin",
    "SM00220": "PKin",
    "PF07714": "PKin",
    "SM00219": "PKin",
    "cd05098": "PKin",
    "cd05078": "PKin",
    "cd14057": "PKin",
    "SM00133": "PKin",
    "PS51285": "PKin",
    "PF00433": "PKin",
    "cd14066": "PKin",
    "PF00085": "TRX",
    "PS51352": "TRX",
    "PF01323": "TRX",
    "cd02947": "TRX",
    "SSF52833": "TRX",
    "G3DSA:3.40.30.10": "TRX",
    "PF03081": "Exo70",
    "PS51514": "BRX",
    "PF08381": "BRX",
    "PF02892": "BED",
    "SM00614": "BED",
    "PS50808": "BED",
    "G3DSA:2.20.25.80": "WRKY",
    "SSF118290": "WRKY",
    "SM00774": "WRKY",
    "PF03106": "WRKY",
    "PS50811": "WRKY",
    "PF00403": "HMA",
    "PS50846": "HMA",
    "cd00371": "HMA",
    "G3DSA:3.30.70.100": "HMA",
    "SSF55008": "HMA",
    "G3DSA:2.100.10.30": "JAC",
    "SSF51101": "JAC",
    "PF16458": "JAC",
    "PF01419": "JAC",
    "PS51752": "JAC",
    "SM00915": "JAC",
    "cd09272": "TRANSPOSON",
    "PF03108": "TRANSPOSON",
    "PF10551": "TRANSPOSON",
    "PF13456": "TRANSPOSON",
    "PF00075": "TRANSPOSON",
    "PS50879": "TRANSPOSON",
    "PF00665": "TRANSPOSON",
    "PF13683": "TRANSPOSON",
    "PF13333": "TRANSPOSON",
    "PS50994": "TRANSPOSON",
    "PF03732": "TRANSPOSON",
    "PF13976": "TRANSPOSON",
    "PF14223": "TRANSPOSON",
    "PF14244": "TRANSPOSON",
    "PF07727": "TRANSPOSON",
    "PF00078": "TRANSPOSON",
    "PS50878": "TRANSPOSON",
    "PF08284": "TRANSPOSON",
    "PF17919": "TRANSPOSON",
    "PF13966": "TRANSPOSON",
    "PF13359": "TRANSPOSON",
    "PF13963": "TRANSPOSON",
    "PF13952": "TRANSPOSON",
    "PF02992": "TRANSPOSON",
    "PF10536": "TRANSPOSON",
    "PF17921": "TRANSPOSON",
    "PF05699": "TRANSPOSON",
    "PF14372": "TRANSPOSON",
    "PF03017": "TRANSPOSON",
    "PF17917": "TRANSPOSON",
    "PF03004": "TRANSPOSON",
    "PF04827": "TRANSPOSON",
    "PF13975": "TRANSPOSON",
    "PF03078": "TRANSPOSON",
    "PF14214": "TRANSPOSON",
    "PF13961": "TRANSPOSON",
    "PF04937": "TRANSPOSON",
    "G3DSA:3.10.10.10": "TRANSPOSON",
    "G3DSA:2.40.70.10": "TRANSPOSON",
    "G3DSA:3.10.20.370": "TRANSPOSON",
    "G3DSA:3.30.70.270": "TRANSPOSON",
    "G3DSA:1.10.340.70": "TRANSPOSON",
    "G3DSA:3.40.395.10": "TRANSPOSON",
    "PS51697": "OTHER",
    "cd02989": "OTHER",
    "SUPERFAMILY": "OTHER",
    "Pfam": "OTHER",
    "Gene3D": "OTHER",
}


class Annotation:
    name = ""
    classification = ""
    start = 0
    end = 0
    evalue = 0.0

    def __init__(self, name, classification, start, end, evalue):
        self.name = name
        self.classification = classification
        self.start = start
        self.end = end
        self.evalue = evalue

    def __lt__(self, other):
        return self.start < other.start


class SequenceAnnotation:
    name = ""
    annotations = {}

    def __init__(self, name):
        self.name = name
        self.annotations = {
            "CC": [],
            "TIR": [],
            "NBARC": [],
            "LRR": [],
            # implement rest in future
        }

    def append(self, annotation):
        if annotation.classification in self.annotations.keys():
            self.annotations[annotation.classification].append(annotation)
        else:
            # do nothing for now
            pass


def merge_and_sort(annotations):
    """
    Merge overlapping annotations for each classification and return a sorted list of annotations.
    """
    merged_annotations = {}
    for classification in annotations:
        merged_annotations[classification] = []
        # get annotations for this classification
        classification_annotations = annotations[classification]
        if len(classification_annotations) == 0:
            continue
        # sort classification annotations by start position
        classification_annotations.sort(key=lambda x: x.start)
        # add first annotation to merged annotations
        merged_annotations[classification].append(classification_annotations[0])
        # merge annotation if it overlaps with the last merged annotation
        # otherwise add it to the list of merged annotations
        for annotation in classification_annotations[1:]:
            if annotation.start <= merged_annotations[classification][-1].end:
                if annotation.end > merged_annotations[classification][-1].end:
                    merged_annotations[classification][-1].end = annotation.end
            else:
                merged_annotations[classification].append(annotation)

    # collapse merged annotations into a single list
    collapsed_annotations = []
    for classification in merged_annotations:
        collapsed_annotations += merged_annotations[classification]

    # sort collapsed annotations by start position
    collapsed_annotations.sort(key=lambda x: x.start)
    return collapsed_annotations


def annotation_string(annotations):
    """
    Return a string representation of a list of annotations.
    """
    annotation_strings = []
    for annotation in annotations:
        annotation_strings.append(f"{annotation.classification}")

    return "-".join(annotation_strings)


def parse_hmmer_table(hmmerfile):
    sequence_annotations = {}
    evalue_threshold = 1e-5
    with open(hmmerfile, "r") as file:
        table_reader = csv.DictReader(file, delimiter="\t", )
        for row in table_reader:
            # in future users may want to use their own annotation files, but work with merged tab-delimited for now

            # need to remove the "." suffix from the annotation name
            # Don't use 6th index, use the 11th! This is the c value, it is the per-domain evalue rather than global

            if row["accession"] not in classifications or float(row["evalue"]) > evalue_threshold:
                continue

            sequence_name = row["sequence"]
            classification = classifications[row["accession"]]

            # not worried about annotation direction for now, so if the end is less than the start, swap them
            if int(row["start"]) < int(row["end"]):
                start = int(row["start"])
                end = int(row["end"])
            else:
                start = int(row["end"])
                end = int(row["start"])

            annotation = Annotation(
                row["accession"], classification, start, end, float(row["evalue"])
            )

            if sequence_name not in sequence_annotations:
                sequence_annotations[sequence_name] = SequenceAnnotation(sequence_name)
            else:
                sequence_annotations[sequence_name].append(annotation)
    return sequence_annotations
