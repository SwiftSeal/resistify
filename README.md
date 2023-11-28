# 🕵️ Resistify

Resistify is a program which classifies plant NLRs by their protein domain and motif architecture.
It is designed to be lightweight - no manual database installations or tricky dependencies here!


## Installation

To get started with Resistify:

`pip install resistify`

Resistify requires `biopython` and `scikit-learn==0.24.2`.
It also requires `hmmsearch` and `jackhmmer` - install these [via conda](https://anaconda.org/bioconda/hmmer) or any other means.
I'd recommend creating an environment for it specifically, as scikit-learn dependencies are a bit busted:

```
mamba create -n resistify python==3.9 pip hmmer
pip install resistify
```

A conda distribution is in progress!

## Usage

To run Resistify:

```
resistify <input.fa> <output_directory>
```

Your `input.fa` should contain the amino acid sequences of your proteins of interest.
Multiline and sequence description fields are allowed.
Stop codons "*" are permitted at the end of sequences - internal stop codons are not.

An `output_directory` will be created which will contain the results of your run:
 - `results.tsv` - A table of the length, classification, and predicted functionality of each sequence, as well as the presence of any MADA motif or CJID domain
 - `motifs.tsv` - A table of all the NLRexpress motifs for each sequence
 - `domains.tsv` - A table of all the domains identified for each sequence
 - `nbarc.fasta` - A fasta file of all the NB-ARC domains identified.

## How does it work?

Resistify is a two step process.

First, all sequences are searched for CC, RPW8, TIR, and NB-ARC domains.
This is used to quickly filter out any non-NLR sequences and identify the primary architecture of each NLR.

Secondly, each potential NLR sequence is scanned for CC, TIR, NB-ARC, and LRR associated motifs via NLRexpress. 
These are used as an additional layer of evidence to reclassify each NLR by predicting LRR domains, and predicting any CC or TIR domains which may have been missed in the initial `hmmsearch`.

Resistify will also search for N-terminal MADA motifs and CJID domains that are common to CNLs and TNLs respectively.

## Contributing

Contributions are greatly appreciated!
If you experience any issues running Resistify, please get in touch via the Issues page.