#!/bin/bash

# Set directory paths
BASE_DIR="resistify_models"
ESM_DIR="$BASE_DIR/esm"
PROSTT5_DIR="$BASE_DIR/prostt5"
PROTT5_DIR="$BASE_DIR/prott5"

# URLs of files to download
ESM_FILES=(
    "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt"
)
PROSTT5_FILES=(
    "https://huggingface.co/Rostlab/ProstT5/resolve/main/config.json"
    "https://huggingface.co/Rostlab/ProstT5/resolve/main/pytorch_model.bin"
    "https://huggingface.co/Rostlab/ProstT5/resolve/main/special_tokens_map.json"
    "https://huggingface.co/Rostlab/ProstT5/resolve/main/spiece.model"
    "https://huggingface.co/Rostlab/ProstT5/resolve/main/tokenizer_config.json"
)
PROTT5_FILES=(
    "https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc/resolve/main/config.json"
    "https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc/resolve/main/pytorch_model.bin"
    "https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc/resolve/main/special_tokens_map.json"
    "https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc/resolve/main/spiece.model"
    "https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc/resolve/main/tokenizer_config.json"
)

# Create directories if they don't exist
mkdir -p "$BASE_DIR" "$ESM_DIR" "$PROSTT5_DIR" "$PROTT5_DIR"

# Function to download files
download_file() {
    local url=$1
    local dest_dir=$2
    local file_name=$(basename "$url")

    # Download if file does not exist
    if [[ ! -f "$dest_dir/$file_name" ]]; then
        echo "Downloading $file_name..."
        curl -L -o "$dest_dir/$file_name" "$url"
    else
        echo "$file_name already exists in $dest_dir. Skipping."
    fi
}

# Download main files to BASE_DIR
for url in "${ESM_FILES[@]}"; do
    download_file "$url" "$ESM_DIR"
done

# Download prostt5 files to PROSTT5_DIR
for url in "${PROSTT5_FILES[@]}"; do
    download_file "$url" "$PROSTT5_DIR"
done

for url in "${PROTT5_FILES[@]}"; do
    download_file "$url" "$PROTT5_DIR"
done

echo "All downloads completed."

