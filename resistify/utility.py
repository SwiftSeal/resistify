import sys
import os
import csv
import logging
import json
import hashlib
import requests
from Bio import SeqIO
from resistify.annotations import Sequence

log = logging.getLogger(__name__)


def create_output_directory(outdir):
    try:
        expanded_outdir = os.path.expanduser(os.path.expandvars(outdir))
        os.makedirs(expanded_outdir, exist_ok=True)
        log.debug(f"Output directory created at {expanded_outdir}")
        return expanded_outdir
    except OSError as e:
        log.error(f"Error creating output directory: {e}")
        sys.exit(1)


def parse_fasta(path):
    sequences = []
    with open(path) as file:
        for record in SeqIO.parse(file, "fasta"):
            # need to remove asterisk, interferes with hmmsearch
            sequence_str = str(record.seq).strip("*")
            if "*" in sequence_str:
                log.error(f"Internal stop codon detected in sequence {record.id}")
                sys.exit(1)
            if len(sequence_str) > 100000:
                log.warning(
                    f"Sequence {record.id} is too long (>100k aa). This sequence will not be analysed!"
                )
                continue
            sequences.append(Sequence(record.id, sequence_str))

    return sequences


def save_fasta(sequences, path, classified_only=False):
    with open(path, "w") as file:
        for sequence in sequences:
            if classified_only and sequence.classification is None:
                continue
            file.write(f">{sequence.id}\n")
            file.write(f"{sequence.seq}\n")
    return path


def result_table(sequences, results_dir, type, retain=False):
    with open(os.path.join(results_dir, "results.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        if type == "nlr":
            table_writer.writerow(
                [
                    "Sequence",
                    "Length",
                    "Motifs",
                    "Domains",
                    "Classification",
                    "NBARC_motifs",
                    "MADA",
                    "MADAL",
                    "CJID",
                ]
            )

            nbarc_motifs = [
                "VG",
                "P-loop",
                "RNSB-A",
                "Walker-B",
                "RNSB-B",
                "RNSB-C",
                "RNSB-D",
                "GLPL",
                "MHD",
            ]

            for sequence in sequences:
                if sequence.type == "NLR" or retain:
                    n_nbarc_motifs = 0
                    for motif in nbarc_motifs:
                        if len(sequence.motifs[motif]) > 0:
                            n_nbarc_motifs += 1

                    table_writer.writerow(
                        [
                            sequence.id,
                            len(sequence.seq),
                            sequence.motif_string(),
                            sequence.domain_string,
                            sequence.classification,
                            n_nbarc_motifs,
                            sequence.mada,
                            sequence.madal,
                            sequence.cjid,
                        ]
                    )
        elif type == "prr":
            table_writer.writerow(
                [
                    "Sequence",
                    "Type",
                    "Classification",
                    "Signal_peptide",
                ]
            )
            for sequence in sequences:
                if sequence.type == "RLP" or sequence.type == "RLK":
                    table_writer.writerow(
                        [
                            sequence.id,
                            sequence.type,
                            sequence.classification,
                            sequence.signal_peptide,
                        ]
                    )


def domain_table(sequences, results_dir):
    with open(os.path.join(results_dir, "domains.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(["Sequence", "Domain", "Start", "End"])
        for sequence in sequences:
            for annotation in sequence.merged_annotations:
                table_writer.writerow(
                    [
                        sequence.id,
                        annotation.domain,
                        annotation.start,
                        annotation.end,
                    ]
                )


def annotation_table(sequences, results_dir):
    with open(os.path.join(results_dir, "annotations.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(
            ["Sequence", "Domain", "Start", "End", "E_value", "Score", "Source"]
        )
        for sequence in sequences:
            for annotation in sequence.annotations:
                table_writer.writerow(
                    [
                        sequence.id,
                        annotation.domain,
                        annotation.start,
                        annotation.end,
                        annotation.evalue,
                        annotation.score,
                        annotation.source,
                    ]
                )


def motif_table(sequences, results_dir):
    from resistify.nlrexpress import MOTIF_SPAN_LENGTHS

    output_path = os.path.join(results_dir, "motifs.tsv")
    with open(output_path, "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(
            [
                "Sequence",
                "Motif",
                "Position",
                "Probability",
                "Downstream_sequence",
                "Motif_sequence",
                "Upstream_sequence",
            ]
        )
        for sequence in sequences:
            for motif in sequence.motifs:
                for item in sequence.motifs[motif]:
                    aa_sequence = sequence.seq
                    downstream_sequence = aa_sequence[item.position - 5 : item.position]
                    motif_sequence = aa_sequence[
                        item.position : item.position + MOTIF_SPAN_LENGTHS[motif]
                    ]
                    upstream_sequence = aa_sequence[
                        item.position + MOTIF_SPAN_LENGTHS[motif] : item.position
                        + MOTIF_SPAN_LENGTHS[motif]
                        + 5
                    ]
                    table_writer.writerow(
                        [
                            sequence.id,
                            motif,
                            item.position,
                            item.probability,
                            downstream_sequence,
                            motif_sequence,
                            upstream_sequence,
                        ]
                    )


def extract_nbarc(sequences, results_dir):
    """
    Extract all nbarc domains of all proteins into a fasta file.
    """
    with open(os.path.join(results_dir, "nbarc.fasta"), "w") as file:
        for sequence in sequences:
            count = 1
            for annotation in sequence.merged_annotations:
                if annotation.domain == "NB-ARC":
                    file.write(f">{sequence.id}_{count}\n")
                    file.write(f"{sequence.seq[annotation.start:annotation.end]}\n")
                    count += 1


def coconat_table(sequences, results_dir):
    output_path = os.path.join(results_dir, "coconat.tsv")
    with open(output_path, "w") as f:
        for sequence in sequences:
            for i, probability in enumerate(sequence.cc_probs):
                f.write(f"{sequence.id}\t{i}\t{probability}\n")

def download_files(base_download_path):
    """
    Download files from a JSON configuration file to specified directories.
    
    Args:
        config_path (str): Path to the JSON configuration file
        base_download_path (str): Base path where model directories will be created
    
    Returns:
        dict: Updated configuration with downloaded file paths and SHA256 sums
    """
    config_path = os.path.join(os.path.dirname(__file__), "data", "model_paths.json")

    # Read the configuration
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # Ensure base download path exists
    os.makedirs(base_download_path, exist_ok=True)
    
    # Process each model's files
    for model_name, model_config in config.items():
        # Create model-specific directory
        model_dir = os.path.join(base_download_path, model_config['directory'])
        os.makedirs(model_dir, exist_ok=True)
        
        # Download each file
        for file_info in model_config['files']:
            url = file_info['url']
            filename = os.path.basename(url)
            local_path = os.path.join(model_dir, filename)
            
            # Download file
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            with open(local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            # Calculate SHA256
            with open(local_path, 'rb') as f:
                file_hash = hashlib.sha256()
                for chunk in iter(lambda: f.read(4096), b""):
                    file_hash.update(chunk)
            
            # Update file info with local path and hash
            file_info['local_path'] = local_path
            file_info['sha256'] = file_hash.hexdigest()
    
    return config

def verify_files(base_download_path: str):
    # Initialize results
    config_path = os.path.join(os.path.dirname(__file__), "data", "model_paths.json")

    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # Check each model's files
    for model_name, model_config in config.items():
        model_dir = os.path.join(base_download_path, model_config['directory'])
        
        # Check if model directory exists
        if not os.path.exists(model_dir):
            log.error(f"Directory not found for {model_name}: {model_dir}")
            sys.exit(1)
        
        # Check each file
        for file_info in model_config['files']:
            url = file_info['url']
            filename = os.path.basename(url)
            local_path = os.path.join(model_dir, filename)
            
            # Check if file exists
            if not os.path.exists(local_path):
                log.error(f"File not found: {local_path}")
                sys.exit(1)
    
    return