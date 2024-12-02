#!/bin/bash
#SBATCH -p gpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --export=ALL

# NLR basic
python -m resistify.main nlr test/some_nlrs.fa -o test/test_out/basic

# NLR Coconat
python -m resistify.main nlr --coconat test/some_nlrs.fa -o test/test_out/coconat

# NLR retain
python -m resistify.main nlr --retain test/uvr8.fa -o test/test_out/retain

# PRR basic
python -m resistify.main prr test/fls2.fa -o test/test_out/prr
