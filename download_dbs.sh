
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.64-96.0/interproscan-5.64-96.0-64-bit.tar.gz -P $TMPDIR
tar -xzf $TMPDIR/interproscan-5.64-96.0-64-bit.tar.gz -C $TMPDIR

cp $TMPDIR/interproscan-5.64-96.0/data/pfam/36.0/* ./resistify/data/pfam/
cp $TMPDIR/interproscan-5.64-96.0/data/smart/9.0/* ./resistify/data/smart/
cp $TMPDIR/interproscan-5.64-96.0/data/superfamily/1.75/* ./resistify/data/superfamily/
cp $TMPDIR/interproscan-5.64-96.0/data/gene3d/4.3.0/* ./resistify/data/gene3d/