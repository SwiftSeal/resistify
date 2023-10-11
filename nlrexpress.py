import subprocess
import sys
from dataclasses import dataclass
import numpy as np
from sklearn.neural_network import MLPClassifier
import pickle
from Bio import SeqIO
import pickle
from sklearn.neural_network import MLPClassifier

motif_span_lengths = {
    "extEDVID": 12,
    "bA": 10,
    "aA": 7,
    "bC": 8,
    "aC": 6,
    "bDaD1": 16,
    "aD3": 13,
    "VG": 5,
    "P-loop": 9,
    "RNSB-A": 10,
    "Walker-B": 8,
    "RNSB-B": 7,
    "RNSB-C": 10,
    "RNSB-D": 9,
    "GLPL": 5,
    "MHD": 3,
    "LxxLxL": 6,
}


@dataclass
class ModelData:
    """
    Base class for prediction models
    """

    model: MLPClassifier
    modelPath: str
    params: dict
    name: str

    def loadModel(self, path: str) -> 'ModelData':
        """ Loads model binary file into the object """
        model = pickle.load(open(path, 'rb'))
        return ModelData(name=self.name, modelPath=path, model=model, params=self.params)

    def train(self, X, Y, params, outFile: str, printToFile=False) -> 'ModelData':
        """ trains model """
        model = MLPClassifier()
        for param in params:
            setattr(model, param, params[param])

        model.fit(X, Y)

        if printToFile:
            pickle.dump(model, open(outFile, 'wb'))

        return ModelData(name=self.name, modelPath=outFile, model=model, params=params)

@dataclass
class ModuleData:
    """
    Base class for prediction modules
    """

    predictors: dict
    modelsPath: dict
    predictorsNames: list

    @classmethod
    def loadModels(cls, modelsPath: dict) -> 'ModuleData':
        """ Loads each model binary file into the object """
        predictors = {}
        predictorNames = []
        for predictorName in modelsPath:
            predictor = ModelData(name=predictorName, model='', modelPath=modelsPath[predictorName], params=[])
            predictors[predictorName] = predictor.loadModel(modelsPath[predictorName])
            predictorNames.append(predictorName)

        return cls(modelsPath=modelsPath, predictors=predictors, predictorsNames=predictorNames)



def parse_jackhmmer(file, iteration = False) -> dict:
    hmm_dict = {}

    header1 = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    with open(file) as f:
        i = 1
        lines = f.readlines()
        while i < len(lines):
            line = lines[i].split()
            if line[0] == "//":
                # if we hit // then we are onto the new entry - skip forward 2 lines to get to next NAME
                i += 2
            elif line[0] == "NAME":
                name = line[1]
                # if this is an iteration, need to strip "-i<iteration number>" from the name
                if iteration:
                    name = name[:-3]
                hmm_dict[name] = []
                # There is an extra field in the iteration file, so we need to skip an extra line
                if iteration == False:
                    i += 14
                else:
                    i += 15
            elif line[0] == "HMM":
                # we have hit the start of the matrix, skip forward 5 lines to get to the matrix
                i += 5
            else:
                # we are in a matrix entry for an amino acid
                # create a dictionary of amino acids to values, then skip forward 3 lines to move to next amino acid
                # ensure that all values are floats
                for value in line[1:21]:
                    try:
                        float(value)
                    except ValueError:
                        print(f"Error: {value} is not a float. {name} {file}")
                        sys.exit(1)
                value_dict = dict(zip(header1, line[1:21]))
                hmm_dict[name].append(value_dict)
                i += 3
    
    return hmm_dict

def generateInputFile( seqData:dict, hmm_it1:dict, hmm_it2:dict) -> dict:
    data = {}

    for name in seqData:
        # get the dna sequence
        seq = seqData[name]
        # make blank list for this sequence
        data[name] = []

        # for each amino acid in the sequence
        for i, aa in enumerate(seq):
            # add a blank list for this amino acid
            data[name].append([])
            for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                data[name][-1].append(val)

            if name in hmm_it2:
                for k, (key, val) in enumerate( hmm_it2[name][i].items() ):
                    data[name][-1].append(val)
            else:
                for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                    data[name][-1].append(val)

    return data

def generate_matrix(sequence_data, motif_span):
    matrix = []
    for sequence in sequence_data:
        sequence_length = len(sequence_data[sequence])
        for i in range(sequence_length):
            # make sure we are within the sequence bounds
            if i >= 5 and i < sequence_length - (motif_span + 5):
                features = input_data[sequence]

                matrix.append([])

                # add the hmm values for this range
                for j in range(-5, motif_span + 6):
                    matrix[-1] += features[i + j]

    matrix = np.array(matrix, dtype=float)
    return matrix


#subprocess.run(["jackhmmer", "--cpu", "6", "-N", "2", "-E", "1e-5", "--domE", "1e-5", "--chkhmm", "test", "sample.fasta", "~/scratch/NLRexpress/hmmer_db/targetDB.fasta"])

jackhmmer_iteration_1 = parse_jackhmmer("test-1.hmm")
jackhmmer_iteration_2 = parse_jackhmmer("test-2.hmm", iteration = True)

sequence_data = {}

with open("sample.fasta") as f:
    # this is replicating the way it's done in nlrexpress code. There may be a better alternative
    for record in SeqIO.parse(f, "fasta"):
        sequence_data[record.id] = record.seq

input_data = generateInputFile(sequence_data, jackhmmer_iteration_1, jackhmmer_iteration_2)

nbs_express = ModuleData.loadModels(
    modelsPath = {
        'VG': 'models/MLP_NBS_VG.pkl',
        'P-loop': 'models/MLP_NBS_P-loop.pkl',
        'RNSB-A': 'models/MLP_NBS_RNSB-A.pkl',
        'RNSB-B': 'models/MLP_NBS_RNSB-B.pkl',
        'RNSB-C': 'models/MLP_NBS_RNSB-C.pkl',
        'RNSB-D': 'models/MLP_NBS_RNSB-D.pkl',
        'Walker-B': 'models/MLP_NBS_Walker-B.pkl',
        'GLPL': 'models/MLP_NBS_GLPL.pkl',
        'MHD': 'models/MLP_NBS_MHD.pkl'
    })

#results = {}

for predictor in nbs_express.predictors:
    motif_span = motif_span_lengths[predictor]
    matrix = generate_matrix(sequence_data, motif_span)
    result = nbs_express.predictors[predictor].model.predict_proba(matrix)
    result_index = 0
    for sequence in sequence_data:
        sequence_length = len(sequence_data[sequence])
        for i in range(len(sequence_data[sequence])):
            # make sure we are within the sequence bounds
            if i >= 5 and i < sequence_length - (motif_span + 5):
                value = round(result[result_index][1], 4)
                if value > 0.8:
                    print(sequence, i, predictor, value)
                result_index += 1

"""
countpos = {predictor: 0 for predictor in results}
for sequence in sequence_data:
    for i in range(len(sequence_data[sequence])):
        for predictor in results:
            if i >= 5 and i < sequence_length - (motif_span + 5):
                value = round(results[predictor][countpos[predictor]][1], 4)
                if value > 0.8:
                    print(sequence, i, predictor, value)
                countpos[predictor] += 1
"""
