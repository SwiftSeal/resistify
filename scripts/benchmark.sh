#!/bin/bash

#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -p gpu
#SBATCH --gpus a100:1
#SBATCH -o benchmark_results/benchmark.out
#SBATCH -e benchmark_results/benchmark.out

set -euo pipefail

# set cache to temporary directory to ensure fresh downloads are pulled
export XDG_CACHE_HOME="$TMPDIR/.cache"

FASTA="$SCRATCH/dm.fa"  # hardcode path here

uv run resistify download # make sure cache aint an issue
export HF_HUB_OFFLINE=1

uv run resistify nlr "$FASTA" -o benchmark_results/nlr
uv run resistify nlr --retain "$FASTA" -o benchmark_results/nlr_retain
uv run resistify nlr --coconat "$FASTA" -o benchmark_results/nlr_coconat
uv run resistify prr "$FASTA" -o benchmark_results/prr
