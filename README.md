# Resistify

Resistify is a lightweight program designed to classify NLRs by their protein domain architecture.
I have created this program as an alternative to several similar programmes for a couple of reasons.
 
The first is to move away from using InterProScan as a dependency.
While InterProScan is a useful resource for annotating protein domains, it very feature-rich and can be challenging to set up on a new system.
It's distribution isn't well supported by conda which is an additional challenge when integrating it into automated workflows.
Resistify comes packaged with all the necessary databases so you don't have to worry about setting them up manually

Secondly, I've created this to be as free of dependencies as possible.
This allows Resistify to be easily distributed and quickly installed!

I'm grateful to the authors of NLRexpress for the motif models used in this program.

## Installation

To get started with Resistify:

`pip install resistify`

Resistify requires `biopython` and `scikit-learn==0.24.2`.
I'll work on proper dependency management soon, but `scikit-learn` makes things pretty tricky - I'd recommend making an environment with `0.24.2` pre-installed rather than letting pip try to handle it itself!
It also requires `hmmsearch` and `jackhmmer` - install these via conda or any other means.
A conda distribution is in progress!

## Usage

To run Resistify:

```
resistify <input.fa> <output_directory>
```

Your `input.fa` should contain the amino acid sequences of your proteins of interest.
Multiline and sequence description fields are allowed.

An `output_directory` will be created which will contain the results of your run:
 - `results.tsv` - A table of the length, classification, and predicted functionality of each sequence
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

## Future improvements

Once the core functionality is stable, I will begin integrating NLR-associated into the pipeline.