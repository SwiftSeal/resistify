<div align="center">

# Resistify 🍃

[![Anaconda-Server Badge](https://anaconda.org/bioconda/resistify/badges/version.svg)](https://anaconda.org/bioconda/resistify)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/resistify/badges/latest_release_date.svg)](https://anaconda.org/bioconda/resistify)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/resistify/badges/downloads.svg)](https://anaconda.org/bioconda/resistify)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19187488.svg)](https://doi.org/10.5281/zenodo.19187488)

</div>

Resistify is a program which rapidly identifies and classifies plant resistance genes from protein sequences.
It is designed to be lightweight and easy to use.

![A screenshot of the help interface of resistify](assets/terminal.png)

## Installation

### Conda

`Resistify` is available via the [Bioconda](https://anaconda.org/bioconda/resistify) channel:

```
conda create -n resistify resistify
conda activate resistify
```

> [!NOTE]
> If you want to use the GPU-accelerated pipelines, conda may fail to install a GPU-ready version of `pytorch`. If this occurs, try adding `pytorch-gpu` to the conda environment.

### Singularity

Containers are also available through the [biocontainers repository](https://quay.io/repository/biocontainers/resistify?tab=tags).
To use these with `singularity`, simply run:

`singularity exec docker://quay.io/biocontainers/resistify:<tag-goes-here> resistify`

## Usage

### Identifying NLRs

To predict NLRs within a set of protein sequences, simply run:

```bash
resistify nlr $PROTEIN_FASTA -o $RESULTS_DIR
```

and `Resistify` will identify and classify NLRs, and return some files:
 - `results.tsv` - A table containing the primary results of `Resistify`.
 - `motifs.tsv` - A table of all the NLR motifs identified for each sequence.
 - `domains.tsv` - A table of all the domains identified for each sequence.
 - `annotations.tsv` - A table of the raw annotations for each sequence.
 - `nbarc.fa` - A fasta file of all the NB-ARC domains identified.
 - `nlr.fa` - A fasta file of all NLRs identified.

By default, `Resistify` will only return sequences with identifiable NB-ARC domains.
If you wish to identify highly fragmented NLRs, you can use the `--retain` option which will predict and report NLR-associated motifs for all sequences.

If you want to increase the sensitivity of coiled-coil domain annotation, you can use the option `--coconat`.
This will use [CoCoNat](https://doi.org/10.1093/bioinformatics/btad495) to predict coiled-coil domains.
In practice, I wouldn't expect this mode to pick up on a significant number of missed CC domains, but it can pick up on cryptic CCs that do not have an identifiable EDVID motif.

#### How does it work?

`Resistify` carries out an initial search for common NLR domains to quickly filter and annotate the input sequences.
Then, `Resistify` executes a re-implementation of `NLRexpress` to conduct a fast and accurate search for NLR-associated motifs.
If `--coconat` is used, this will also be executed to scavenge for potentially missed coiled-coil domains.
Together, this evidence is used to classify NLRs according to their domain architecture.

### Identifying PRRs

> [!IMPORTANT]
> This pipeline is currently in development - due to other commitments I can't currently benchmark this properly and make no guarantees to its accuracy yet! Feedback is appreciated.

To predict PRRs within a set of protein sequences, simply run:

```
resistify prr $PROTEIN_FASTA -o $RESULTS_DIR
```

and `Resistify` will identify and classify PRRs, and return some files:
 - `results.tsv` - A table containing the primary results of `Resistify`.
 - `motifs.tsv` - A table of all the LRR motifs identified for each sequence.
 - `domains.tsv` - A table of all the domains identified for each sequence.
 - `annotations.tsv` - A table of the raw annotations for each sequence.
 - `prr.fa` - A fasta file of all PRRs identified.

> [!WARNING]
> This pipeline is GPU-accelerated and will be slow on CPU only.

#### How does it work?

First, `Resistify` searches for domains associated with a [recently described classification system](https://doi.org/10.1016/j.molp.2024.02.014) for RLP/RLKs.
Then, a re-implementation of [`TMbed`](https://github.com/BernhoferM/TMbed) is used to predict transmembrane domains - sequences with a single α-helix transmembrane domain and an extracellular domain of at least 50 amino acids are considered as RLPs.
Finally, `NLRexpress` is used to identify LRR domains.

Sequences are classified as being either RLPs or RLKs depending on the presence of an internal kinase domain, and are classified according to their extracellular domain.

### Downloading model data

Models are downloaded automatically at runtime.
If you're going to run `Resistify` frequently, for example as part of a pipeline, you might want to pre-download models first.
This can be done via:

```bash
resistify download
```

This downloads models to the default Hugging Face / PyTorch Hub cache (`$HOME/.cache`, or `$XDG_CACHE_HOME` if set).
You can then set `$HF_HUB_OFFLINE=1` which will speed up `Resistify` slightly and prevents too many requests being sent:

```bash
HF_HUB_OFFLINE=1 resistify nlr $PROTEIN_FASTA -o $RESULTS_DIR
```

## Results

### results.tsv (nlr)

| Sequence | Length | LRR_Length | Motifs                  | Domains | Classification | NBARC_motifs | MADA  | CJID  |
| -------- | -------| ---------- | ----------------------- | ------- | -------------- | ------------ | ----- | ----- |
| ZAR1     | 852    | 306        | CNNNNNNNNNLLLLLLLLLLLLL | mCNL    | CNL            | 9            | False | False |

The main column of interest is "Classification", where we can see that it has been identified as a canonical CNL.
The "Motifs" column indicates the series of NLR-associated motifs identified across the sequence - this can be useful if an NLR has an undetermined or unexpected classification.
The columns "MADA", and "CJID" correspond to common NLR sequence signatures.

### results.tsv (prr)

| Sequence | Length | Extracellular_Length | LRR_Length | Type | Classification | Signal_peptide |
| --- | --- | --- | --- | --- | --- | --- |
| fls2 | 1173 | 806 | 675 | RLK | LRR | True |

For PRRs, sequences can be of the type RLP or RLK - both are single pass transmembrane proteins, and RLKs have an internal kinase domain.
Classification refers to the domains identified in the external region.
If multiple domains are identified, they will each be reported as a semi-colon separated list.
If a signal peptide is identified in the sequence, this is reported accordingly.

### motifs.tsv

| Sequence | Motif    | Position | Probability | Downstream_sequence | Motif_sequence | Upstream_sequence |
|----------|----------|----------|-------------|---------------------|----------------|-------------------|
| ZAR1     | extEDVID | 66       | 0.85        | LVADL               | RELVYEAEDILV   | DCQLA             |
| ZAR1     | VG       | 160      | 1.0         | YDHTQ               | VVGLE          | GDKRK             |
| ZAR1     | P-loop   | 189      | 1.0         | IMAFV               | GMGGLGKTT      | IAQEV             |
| ZAR1     | RNBS-A   | 212      | 1.0         | EIEHR               | FERRIWVSVS     | QTFTE             |
| ZAR1     | Walker-B | 260      | 1.0         | QYLLG               | KRYLIVMD       | DVWDK             |
| ZAR1     | RNBS-B   | 291      | 1.0         | RGQGG               | SVIVTTR        | SESVA             |
| ZAR1     | RNBS-C   | 318      | 1.0         | HRPEL               | LSPDNSWLLF     | CNVAF             |
| ZAR1     | GLPL     | 357      | 1.0         | VTKCK               | GLPLT          | IKAVG             |
| ZAR1     | RNBS-D   | 418      | 1.0         | SHLKS               | CILTLSLYP      | EDCVI             |
| ZAR1     | MHD      | 487      | 1.0         | IITCK               | IHD            | MVRDL             |
| ZAR1     | LxxLxL   | 512      | 0.96        | PEGLN               | CRHLGI         | SGNFD             |
| ZAR1     | LxxLxL   | 532      | 0.89        | KVNHK               | LRGVVS         | TTKTG             |
| ZAR1     | LxxLxL   | 561      | 1.0         | TDCKY               | LRVLDI         | SKSIF             |
| ZAR1     | LxxLxL   | 588      | 1.0         | ASLQH               | LACLSL         | SNTHP             |
| ZAR1     | LxxLxL   | 612      | 1.0         | EDLHN               | LQILDA         | SYCQN             |
| ZAR1     | LxxLxL   | 636      | 1.0         | VLFKK               | LLVLDM         | TNCGS             |
| ZAR1     | LxxLxL   | 660      | 1.0         | GSLVK               | LEVLLG         | FKPAR             |
| ZAR1     | LxxLxL   | 686      | 1.0         | KNLTN               | LRKLGL         | SLTRG             |
| ZAR1     | LxxLxL   | 713      | 1.0         | INLSK               | LMSISI         | NCYDS             |
| ZAR1     | LxxLxL   | 741      | 1.0         | TPPHQ               | LHELSL         | QFYPG             |
| ZAR1     | LxxLxL   | 766      | 1.0         | HKLPM               | LRYMSI         | CSGNL             |
| ZAR1     | LxxLxL   | 793      | 1.0         | NTHWR               | IEGLML         | SSLSD             |
| ZAR1     | LxxLxL   | 818      | 1.0         | QSMPY               | LRTVTA         | NWCPE             |

Here, the positions, probabilities, and sequence of NLRexpress motif hits are listed.
The five amino acids upstream and downstream of the motif site are also provided.
In PRR mode, only LRR motifs will be reported.

### domains.tsv

| Sequence | Domain | Start | End |
|----------|--------|-------|-----|
| ZAR1     | MADA   | 1     | 21  |
| ZAR1     | CC     | 5     | 128 |
| ZAR1     | NBARC  | 163   | 337 |
| ZAR1     | LRR    | 512   | 817 |

This file contains the coordinates of the domains identified by `Resistify`.

### annotations.tsv

| Sequence | Domain | Start | End | Accession  | Score  | Source    |
|----------|--------|-------|-----|------------|--------|-----------|
| ZAR1     | MADA   | 1     | 21  | CC_motif_1 | 16.2   | hmmer     |
| ZAR1     | CC     | 5     | 128 | PF18052    | 69.72  | hmmer     |
| ZAR1     | NBARC  | 163   | 337 | PF00931    | 196.02 | hmmer     |
| ZAR1     | LRR    | 512   | 817 |            |        | resistify |
| ZAR1     | LRR    | 541   | 781 | PF23598    | 99.53  | hmmer     |
| ZAR1     | LRR    | 676   | 808 | PF25019    | 31.39  | hmmer     |

This file contains the raw annotations for each sequence, and the method which was used to identify them.

### plots.tar.gz

By default, `Resistify` generates some rudimentary plots for each protein, saved as a `plots.tar.gz` archive.
You can disable these via `--no-draw` if ya want.

![An SVG of ZAR1](assets/zar1.svg)

## Frequently asked questions

**Q: Can `Resistify` be used to predict resistance genes from genomic data?**

**A:** Unfortunately, `Resistify` cannot be directly applied to a genome to predict resistance genes, unlike tools such as `NLR-Annotator`.
If gene annotations are unavailable for your genome, my advice would be to use a tool like [`Helixer`](https://github.com/weberlab-hhu/Helixer) or [`ANNEVO`](https://github.com/xjtu-omics/ANNEVO) to perform *ab initio* gene prediction first, then pass these to `Resistify`.
Currently, I find that `Helixer` tends to identify more NLRs than `ANNEVO` (in *Solanum*):

![A barplot of the number of NLRs identified by Helixer vs ANNEVO](assets/predictions.png)

**Q: According to the Motif string, some of my genes have NLR motifs in unexpected places - are these significant?**

**A:** False positives do occur for the motif predictions, and unexpected predictions such as a single CC motif in the LRR domain are unlikely to be representative of a true domain annotation.
False positives shouldn't interfere with the classification accuracy.

## Contributing

Contributions are greatly appreciated!
If you experience any issues running Resistify, please get in touch via the Issues page.
If you have any suggestions for additional features, get in touch!

## Citing

> Smith M., Jones J. T., Hein I. (2025) Resistify: A Novel NLR Classifier That Reveals Helitron-Associated NLR Expansion in Solanaceae. *Bioinformatics and Biology Insights*. 2025;19. [doi:10.1177/11779322241308944](https://doi.org/10.1177/11779322241308944)

You must also cite:

> Martin, E. C., Spiridon, L., Goverse, A., & Petrescu, A. J. (2022). NLRexpress—A bundle of machine learning motif predictors—Reveals motif stability underlying plant Nod-like receptors diversity. *Frontiers in Plant Science*, 13, 975888. https://doi.org/10.3389/fpls.2022.975888

If you use the `CoCoNat` module, please cite:

> Madeo, G., Savojardo, C., Manfredi, M., Martelli, P. L., & Casadio, R. (2023). CoCoNat: a novel method based on deep learning for coiled-coil prediction. *Bioinformatics*, 39(8), btad495. https://doi.org/10.1093/bioinformatics/btad495

If you use the PRR module, please cite:

> Bernhofer, M., & Rost, B. (2022). TMbed: transmembrane proteins predicted through language model embeddings. *BMC bioinformatics*, 23(1), 326. https://doi.org/10.1186/s12859-022-04873-x
