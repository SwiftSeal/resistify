from resistify._loguru import logger
from resistify.annotations import Sequence, NBARC_MOTIFS
import os
import csv

def collect_data(queries, results_dir):
    sequences = {}

    logger.debug("Collecting data from results directory...")
    with open(os.path.join(results_dir, "results.tsv")) as results_file:
        reader = csv.DictReader(results_file, delimiter="\t")
        for row in reader:
            if row["Sequence"] in queries:
                sequences[row["Sequence"]] = {
                    "length": int(row["Length"]),
                    "domains": [],
                    "motifs": []
                }

    with open(os.path.join(results_dir, "domains.tsv")) as domains_file:
        reader = csv.DictReader(domains_file, delimiter="\t")
        for row in reader:
            if row["Sequence"] in sequences:
                sequences[row["Sequence"]]["domains"].append({
                    "start": int(row["Start"]),
                    "end": int(row["End"]),
                    "domain": row["Domain"]
                })
    
    with open(os.path.join(results_dir, "motifs.tsv")) as motifs_file:
        reader = csv.DictReader(motifs_file, delimiter="\t")
        for row in reader:
            if row["Sequence"] in sequences:
                sequences[row["Sequence"]]["motifs"].append({
                    "position": int(row["Position"]),
                    "motif": row["Motif"]
                })
    
    return sequences

def draw(args):
    queries = [q.strip() for q in args.query.split(",") if q.strip()]
    sequences = collect_data(queries, args.results)
    # For now, just print the collected data for verification
    logger.debug(sequences)
