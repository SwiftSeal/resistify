#!/bin/bash

#SBATCH -c 16
#SBATCH --mem 32G
#SBATCH -p gpu
#SBATCH --gpus a100:1
#SBATCH -o benchmark_results/benchmark.out
#SBATCH -e benchmark_results/benchmark.out

set -euo pipefail

FASTA="$SCRATCH/dm.pep.fa"  # hardcode path here

uv run resistify nlr "$FASTA" -o benchmark_results/nlr
uv run resistify nlr --retain "$FASTA" -o benchmark_results/nlr_retain
uv run resistify nlr --coconat "$FASTA" -o benchmark_results/nlr_coconat
uv run resistify prr "$FASTA" -o benchmark_results/prr
