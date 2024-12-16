#!/bin/bash
#SBATCH -J resistify_nlr_test
#SBATCH -p long
#SBATCH -c 16
#SBATCH --mem=16G

source activate resistify_dev

curl -o $TMPDIR/refplantnlr.fa -LO "https://doi.org/10.1371/journal.pbio.3001124.s013"
curl -o $TMPDIR/dm.fa.gz -O "https://spuddb.uga.edu/data/DM_1-3_516_R44_potato.v6.1.hc_gene_models.pep.fa.gz"

python -m resistify.main nlr $TMPDIR/refplantnlr.fa -o test_nlr_refplantnlr
python -m resistify.main nlr $TMPDIR/dm.fa.gz -o test_nlr_dm

python -m resistify.main nlr --coconat $TMPDIR/refplantnlr.fa -o test_nlr_refplantnlr_coconat
python -m resistify.main nlr --coconat $TMPDIR/dm.fa.gz -o test_nlr_dm_coconat