#!/bin/bash
#SBATCH -J resistify_prr_test
#SBATCH -p gpu
#SBATCH -c 16
#SBATCH --mem=16G

curl -o dm.fa.gz -O "https://spuddb.uga.edu/data/DM_1-3_516_R44_potato.v6.1.hc_gene_models.pep.fa.gz"

python -m resistify.main prr dm.fa.gz -o test_nlr_dm