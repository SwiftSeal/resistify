wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.64-96.0/interproscan-5.64-96.0-64-bit.tar.gz -P $TMPDIR
tar -xzf $TMPDIR/interproscan-5.64-96.0-64-bit.tar.gz -C $TMPDIR

cp $TMPDIR/interproscan-5.64-96.0/data/pfam/36.0/pfam_a.hmm ./resistify/data/pfam.hmm
cp $TMPDIR/interproscan-5.64-96.0/data/smart/9.0/smart.HMMs ./resistify/data/smart.hmm
cp $TMPDIR/interproscan-5.64-96.0/data/superfamily/1.75/hmmlib_1.75 ./resistify/data/superfamily.hmm
cp $TMPDIR/interproscan-5.64-96.0/data/gene3d/4.3.0/gene3d_main.hmm ./resistify/data/gene3d.hmm
cp $TMPDIR/interproscan-5.64-96.0/data/gene3d/4.3.0/model_to_family_map.tsv ./resistify/data/gene3d.tsv

rm -rf $TMPDIR/interproscan-5.64-96.0-64-bit.tar.gz $TMPDIR/interproscan-5.64-96.0