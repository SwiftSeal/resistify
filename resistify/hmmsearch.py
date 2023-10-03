import subprocess
import os
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

    # superfamily and gene3d don't have properly formatted accession names.
    # So if the source is either of those, we need to reformat these.
    if source == "superfamily":
        logging.info("😊 Fixing superfamily accession names...")
        print_superfamily_accessions(tmp_file)
    elif source == "gene3d":
        logging.info("😊 Fixing gene3d accession names...")
        print_gene3d_accessions(tmp_file)
    elif source == "pfam":
        logging.info("😊 Fixing pfam accession names...")
        print_pfam_accessions(tmp_file)
    elif source == "smart":
        # Don't need to do anything to smart hmmsearch output, just remove whitespace and headers
        with open(tmp_file, "r") as file:
            for line in file:
                if not line.startswith("#"):
                    line = line.split()
                    print("\t".join(line))
    elif source == "cjid":
        # Rename accession to "CJID"
        with open(tmp_file, "r") as file:
            for line in file:
                if not line.startswith("#"):
                    line = line.split()
                    line[4] = "CJID"
                    print("\t".join(line))


def print_superfamily_accessions(superfamily_file):
    """
    Superfamily accession names are not formatted correctly for resistify.
    This function fixes the accession names and prints them to stdout.
    """
    with open(superfamily_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                line = line.split()
                # Accession not provided, but just append SSF
                line[4] = "SSF" + line[3]
                print("\t".join(line))


def print_gene3d_accessions(
    gene3d_file, model_to_family_map="./resistify/data/gene3d.tsv"
):
    """
    Gene3D accession names are not formatted correctly for resistify.
    This function fixes the accession names and prints them to stdout.
    To do this, the model_to_family_map.tsv file provided by interproscan is used.
    """
    model_to_family_map_dict = {}
    with open(model_to_family_map, "r") as file:
        for line in file:
            line = line.split("\t")
            model_to_family_map_dict[line[0]] = line[1]

    with open(gene3d_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                line = line.split()
                key = line[3].split("-")[0]
                if key not in model_to_family_map_dict:
                    logging.warning(f"🤔 {key} not found in gene3d model to family map!")
                    continue
                line[4] = "G3DSA:" + model_to_family_map_dict[key]
                print("\t".join(line))


def print_pfam_accessions(pfam_file):
    with open(pfam_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                line = line.split()
                line[4] = line[4].split(".")[0]
                print("\t".join(line))


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
