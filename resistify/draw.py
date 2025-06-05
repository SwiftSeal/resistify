from resistify._loguru import logger
import os
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches

COLOUR_PALETTE = {
    "CC": "#648FFF",
    "RPW8": "#785EF0",
    "TIR": "#DC267F",
    "NB-ARC": "#FE6100",
    "LRR": "#FFB000",
    "signal_peptide": "#648FFF",
    "alpha_inwards": "#DC267F",
    "PKinase": "#FE6100",
}

MOTIF_IDENTIFIERS = {
    'aA': 'TIR',
    'aC': 'TIR',
    'aD3': 'TIR',
    'bA': 'TIR',
    'bC': 'TIR',
    'bDaD1': 'TIR',
    'extEDVID': 'CC',
    'VG': 'NB-ARC',
    'P-loop': 'NB-ARC',
    'RNSB-A': 'NB-ARC',
    'Walker-B': 'NB-ARC',
    'RNSB-B': 'NB-ARC',
    'RNSB-C': 'NB-ARC',
    'RNSB-D': 'NB-ARC',
    'GLPL': 'NB-ARC',
    'MHD': 'NB-ARC',
    'LxxLxL': 'LRR',
}

def collect_data(queries, results_dir):
    sequences = {}

    logger.debug("Collecting data from results directory...")
    with open(os.path.join(results_dir, "results.tsv")) as results_file:
        reader = csv.DictReader(results_file, delimiter="\t")
        for row in reader:
            # If queries is None, collect all sequences; otherwise, filter by queries
            if queries is None or row["Sequence"] in queries:
                sequences[row["Sequence"]] = {
                    "length": int(row["Length"]),
                    "domains": [],
                    "motifs": []
                }

    with open(os.path.join(results_dir, "domains.tsv")) as domains_file:
        reader = csv.DictReader(domains_file, delimiter="\t")
        for row in reader:
            if row["Sequence"] in sequences:
                #special case - treat MADA as a motif
                if row["Domain"] == "MADA":
                    sequences[row["Sequence"]]["motifs"].append({
                        "position": int(row["Start"]),
                        "motif": "MADA"
                    })
                else:
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

def draw_nlr(args, sequence_data):
    logger.info(f"Drawing NLRs to {args.output}...")
    fig, ax = plt.subplots(figsize=(args.width, args.height), dpi = 600)

    # Calculate y-offsets for multiple sequences to stack them
    y_offset = 0
    y_tick_positions = []
    y_tick_labels = []

    for seq_id, features in sequence_data.items():
        length = features['length']
        domains = features.get('domains', [])
        motifs = features.get('motifs', [])

        # Draw the horizontal line for sequence length
        ax.hlines(y_offset, 0, length, color='black', linestyle='-', linewidth=2, label='Sequence Length' if y_offset == 0 else "", zorder = 1)

        # Draw lollipop for motifs
        if not args.hide_motifs:
            for motif_info in motifs:
                position = motif_info['position']
                motif_name = motif_info['motif']
                motif_type = MOTIF_IDENTIFIERS.get(motif_name)
                motif_color = COLOUR_PALETTE.get(motif_type, 'gray')
                # Draw vertical line ("lollipop stick")
                ax.vlines(position, y_offset, y_offset + 0.15, color='black', linewidth=1, label='Motif' if y_offset == 0 else "", zorder = 2)
                # Place motif name so it starts at the end of the lollipop stick
                ax.text(
                    position - 5,
                    y_offset + 0.15,  # exactly at the end of the vline
                    motif_name,
                    rotation=45,
                    ha='left',       # start from the vline
                    va='bottom',     # align bottom to the vline end
                    fontsize=7,
                    color=motif_color
                )

        # Draw rectangles for domains
        for domain_info in domains:
            start = domain_info['start']
            end = domain_info['end']
            domain_name = domain_info['domain']
            domain_color = COLOUR_PALETTE.get(domain_name, 'gray')
            rect = patches.Rectangle((start, y_offset - 0.1), end - start, 0.2,
                                     facecolor=domain_color, edgecolor='black', zorder = 3)
            ax.add_patch(rect)
            # Add domain name centered in the rectangle
            ax.text(
                (start + end) / 2,
                y_offset,
                domain_name,
                ha='center',
                va='center',
                fontsize=8,
                color='black',
                zorder=4
            )

        y_tick_positions.append(y_offset)
        y_tick_labels.append(seq_id)
        y_offset += 1  # Increment offset for the next sequence

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels(y_tick_labels)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('Position (aa)')
    ax.set_xlim(left=0)
    ax.set_ylim(-0.5, y_offset - 0.5)
    plt.tight_layout()
    plt.savefig(args.output, transparent=True)

def draw_prr(args, sequence_data):
    logger.info(f"Drawing PRRs to {args.output}...")
    fig, ax = plt.subplots(figsize=(args.width, args.height), dpi = 600)

    # Calculate y-offsets for multiple sequences to stack them
    y_offset = 0
    y_tick_positions = []
    y_tick_labels = []

    for seq_id, features in sequence_data.items():
        length = features['length']
        domains = features.get('domains', [])
        motifs = features.get('motifs', [])

        # Draw the horizontal line for sequence length
        # For PRRs, offset the x axis to the start of the transmembrane domain
        for domain in domains:
            if domain['domain'] == 'alpha_inwards':
                x_offset = domain['start'] + (domain['end'] - domain['start'])/2
                logger.debug(f"Transmembrane domain found at {x_offset} for sequence {seq_id}.")
                break
        
        ax.hlines(
            y_offset,
            xmin = x_offset,
            xmax = x_offset - length,
            color='black',
            linestyle='-',
            linewidth=2,
        )

        if not args.hide_motifs:
            for motif_info in motifs:
                position = motif_info['position']
                motif_name = motif_info['motif']
                motif_type = MOTIF_IDENTIFIERS.get(motif_name)
                motif_color = COLOUR_PALETTE.get(motif_type, 'gray')
                # Draw vertical line ("lollipop stick")
                ax.vlines(x_offset - position, y_offset, y_offset + 0.15, color='black', linewidth=1, label='Motif' if y_offset == 0 else "", zorder = 2)
                # Place motif name so it starts at the end of the lollipop stick
                ax.text(
                    x_offset - position - 5,
                    y_offset + 0.15,  # exactly at the end of the vline
                    motif_name,
                    rotation=45,
                    ha='left',       # start from the vline
                    va='bottom',     # align bottom to the vline end
                    fontsize=7,
                    color=motif_color
                )

        # Draw rectangles for domains
        for domain_info in domains:
            start = x_offset - domain_info['start']
            end = x_offset - domain_info['end']
            domain_name = domain_info['domain']
            logger.debug(f"Drawing domain {domain_name} from {start} to {end} for sequence {seq_id}.")
            domain_color = COLOUR_PALETTE.get(domain_name, 'gray')
            rect = patches.Rectangle((min(start, end), y_offset - 0.1), abs(end - start), 0.2,
                                     facecolor=domain_color, edgecolor='black', zorder = 3)
            ax.add_patch(rect)

            # Relabel domain names for PRRs
            if domain_name == "alpha_inwards":
                domain_name = "TM"
            elif domain_name == "signal_peptide":
                domain_name = "SP"
            elif domain_name == "Pkinase":
                domain_name = "Kinase"

            ax.text(
                (start + end) / 2,
                y_offset,
                domain_name,
                ha='center',
                va='center',
                fontsize=8,
                color='black',
                zorder=4
            )

        y_tick_positions.append(y_offset)
        y_tick_labels.append(seq_id)
        y_offset += 1  # Increment offset for the next sequence

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels(y_tick_labels)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('Position (aa)')
    ax.set_ylim(-0.5, y_offset - 0.5)
    plt.tight_layout()
    plt.savefig(args.output, transparent=True)

    

def draw(args):
    if not os.path.exists(args.results_dir):
        logger.error(f"Results directory {args.results_dir} does not exist.")
        return
    
    if args.query is None:
        logger.info("No specific queries provided, plotting all sequences.")
        queries = None
    else:
        queries = [q.strip() for q in args.query.split(",") if q.strip()]

    sequence_data = collect_data(queries, args.results_dir)

    # First, need to detect if we are drawing NLRs or PRRs
    # Can do this by checking for prr.fasta or nlr.fasta in directory
    if "nlr.fasta" in os.listdir(args.results_dir):
        logger.info("Detected NLR sequences, drawing NLRs...")
        draw_nlr(args, sequence_data)
    elif "prr.fasta" in os.listdir(args.results_dir):
        logger.info("Detected PRR sequences, drawing PRRs...")
        draw_prr(args, sequence_data)
