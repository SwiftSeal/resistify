# *Resistify*

*Resistify* is a lightweight and fast program designed to classify NLRs by their protein domain architecture.
It does not require any external databases or manual configuration.

## Installation

To get started with *Resistify*:

`pip install resistify`

*Resistify* requires `biopython` and `scikit-learn==0.24.2`.
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

An `output_directory` will be created which will contain the results of your run:
 - `results.tsv` - A table of the length, classification, and predicted functionality of each sequence, as well as the presence of any MADA motif or CJID domain
 - `motifs.tsv` - A table of all the NLRexpress motifs for each sequence
 - `domains.tsv` - A table of all the domains identified for each sequence
 - `nbarc.fasta` - A fasta file of all the NB-ARC domains identified.

## How does it work?

Resistify is a two step process.

First, all sequences are searched for CC, RPW8, TIR, and NB-ARC domains.
This is used to quickly filter out any non-NLR sequences and identify the primary architecture of each NLR.

Secondly, each potential NLR sequence is scanned for CC, NB-ARC, and LRR associated motifs via the NLRexpress models.
These are used as an additional layer of evidence to reclassify each NLR by predicting LRR domains, and predicting any CC domains which may have been missed in the initial `hmmsearch` which can be less sensitive for this domain.
The functionality of each NLR is predicted by counting the number of conserved NB-ARC motifs.
Currently, any order is accepted (this may change in the future!).

Resistify will also search for N-terminal MADA motifs and CJID domains that are common to CNLs and TNLs respectively.