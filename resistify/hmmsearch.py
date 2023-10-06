import subprocess
import os
import csv
import requests
import tarfile
import shutil
import logging
import tempfile


def hmmsearch(input_fasta, source, database_path, e_value="0.00001", num_cpu="2"):
    # create a temporary file to store the hmmsearch --domtblout output
    with tempfile.NamedTemporaryFile(mode="w") as tmp:
        tmp_file = tmp.name

    cmd = [
        "hmmsearch",
        "--noali",
        "-E",
        e_value,
        "--cpu",
        num_cpu,
        "--domtblout",
        tmp_file,
        database_path,
        input_fasta,
    ]

    logging.debug(f"Running hmmsearch with command: {cmd}")

    logging.info(f"😊 Running hmmsearch for {source}...")
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        logging.info(f"😊 {source} hmmsearch completed successfully...")
    except subprocess.CalledProcessError as e:
        logging.error(f"😞 Error running hmmsearch. Stdout of hmmsearch: {e.stdout}")

    return (source, tmp_file)

def save_fixed_accession(results, results_dir):
    # create temporary file to store the unsorted output
    tmp = tempfile.NamedTemporaryFile(mode="w")
    for result in results:
        logging.info(f"😊 Fixing {result[0]} accession names...")
        model_to_family_map_dict = None
        with open(result[1]) as file:
            for line in file:
                if not line.startswith("#"):
                    line = line.split()
                    if result[0] == "superfamily":
                        # Accession not provided, but just append SSF
                        line[4] = "SSF" + line[3]
                    elif result[0] == "gene3d":
                        # Load model to family map if not already loaded
                        if model_to_family_map_dict is None:
                            model_to_family_map_dict = parse_gene3d_table("./resistify/data/gene3d.tsv")
                        key = line[3].split("-")[0]
                        if key not in model_to_family_map_dict:
                            logging.warning(f"🤔 {key} not found in gene3d model to family map!")
                            continue
                        line[4] = "G3DSA:" + model_to_family_map_dict[key]
                    elif result[0] == "pfam":
                        line[4] = line[4].split(".")[0]
                    # Write the line to the temporary file
                    tmp.write("\t".join(line) + "\n")
    tmp.close()
    # Sort the temporary file and save it to the results directory
    logging.info(f"😊 Saving hmmsearch results...")
    with open(tmp.name, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        sorted_reader = sorted(reader, key=lambda row: row[0])
        with open(os.path.join(results_dir, "hmmsearch_result.tsv"), "w") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerows(sorted_reader)
    
    return os.path.join(results_dir, "hmmsearch_result.tsv")

def parse_gene3d_table(gene3d_file):
    model_to_family_map_dict = {}
    with open(gene3d_file, "r") as file:
        for line in file:
            line = line.split()
            model_to_family_map_dict[line[0]] = line[1]
    return model_to_family_map_dict
        

def get_interproscan_data():
    # Define URLs and paths
    base_url = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.64-96.0/"
    tar_file_url = base_url + "interproscan-5.64-96.0-64-bit.tar.gz"
    tmp_dir = "/tmp"  # You can change this to your desired temporary directory
    interproscan_dir = os.path.join(tmp_dir, "interproscan-5.64-96.0")

    # Download the tar.gz file
    tar_file_path = os.path.join(tmp_dir, "interproscan-5.64-96.0-64-bit.tar.gz")
    response = requests.get(tar_file_url)
    response.raise_for_status()
    with open(tar_file_path, "wb") as f:
        f.write(response.content)

    # Extract the tar.gz file
    with tarfile.open(tar_file_path, "r:gz") as tar:
        tar.extractall(tmp_dir)

    # Copy specific files to the desired location
    data_dir = os.path.join(interproscan_dir, "data")
    resistify_data_dir = "data"

    file_mappings = [
        ("pfam/36.0/pfam_a.hmm", "pfam.hmm"),
        ("smart/9.0/smart.HMM", "smart.hmm"),
        ("superfamily/1.75/hmmlib_1.75", "superfamily.hmm"),
        ("gene3d/4.3.0/gene3d_main.hmm", "gene3d.hmm"),
        ("gene3d/4.3.0/model_to_family_map.tsv", "gene3d.tsv"),
    ]

    for src_filename, dest_filename in file_mappings:
        src_path = os.path.join(data_dir, src_filename)
        dest_path = os.path.join(resistify_data_dir, dest_filename)
        shutil.copy(src_path, dest_path)

    # Clean up temporary files and directories
    os.remove(tar_file_path)
    shutil.rmtree(interproscan_dir)
