import subprocess
from Bio import SeqIO
from time import sleep

def parse_jackhmmer(file) -> dict:
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
                hmm_dict[name] = []
                i += 14
            elif line[0] == "HMM":
                # we have hit the start of the matrix, skip forward 5 lines to get to the matrix
                i += 5
            else:
                # we are in a matrix entry for an amino acid
                # create a dictionary of amino acids to values, then skip forward 3 lines to move to next amino acid
                value_dict = dict(zip(header1, line[1:21]))
                hmm_dict[name].append(value_dict)
                i += 3

def generateInputFile( seqData:dict, hmm_it1:dict, hmm_it2:dict) -> dict:
    data = {}

    for name in seqData:
        seq = seqData[name]
        data[name] = []

        for i, aa in enumerate(seq):
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

#subprocess.run(["jackhmmer", "--cpu", "6" "-E", "1e-5", "--domE", "1e-5", "--chkhmm", "test", "sample.fasta", "~/scratch/NLRexpress/hmmer_db/targetDB.fasta"])

jackhmmer_iteration_1 = parse_jackhmmer("test-1.hmm")

with open("sample.fasta") as f:
    for record in SeqIO.parse(f, "fasta"):
        seq = record.seq