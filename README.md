<img src="https://github.com/RuolinHe/NRPS-motif-Finder/assets/76482251/5ee2d1fc-b026-4be3-bfdc-6ac6d6151a8b" alt="Striking Image" width="600" height="600">

The [cover](https://journals.plos.org/ploscompbiol/issue?id=10.1371/issue.pcbi.v19.i05) of PLOS Computational Biology in 2023 May for this work.

# What is NRPS-motif-Finder
NRPS-motif-Finder is a tool for standardization of Non-ribosomal peptide synthetase (NRPS). It partitions the input NRPS protein sequence by locating these conserved motif, to output a motif-and-intermotif architecture that feeds in subsequent analysis such as C domain classification, NRPS re-engineering.

<img src="https://github.com/RuolinHe/NRPS-motif-Finder-matlab-version/assets/76482251/a0358bd6-3431-4bda-9477-20898debc41a" alt="logo" width="388.9" height="200">

# Supported domains and motifs
Adenylation (A) domain has 12 domain: Aalpha, A1-A5, G-motif, A6-A10. Among them, Aalpha and G-motif were two new motifs proposed in [our paper](https://doi.org/10.1371/journal.pcbi.1011100).

Condensation (C) domain has 10 domain: C1-C10.

Thiolation (T) domain has 2 domain: Talpha, T1. Talpha was one new motif [our paper](https://doi.org/10.1371/journal.pcbi.1011100).

Thioesterase (TE) domain has 1 domain: TE1.

Epimerization (E) domain has 7 domains: E1-E7.

# C domain subtype
One of the most important features is that NRPS-motif-Finder supports the full subtype classification of C domain.

<img src="https://github.com/RuolinHe/NRPS-motif-Finder-matlab-version/assets/76482251/dfc1f54c-b746-4dc1-b9ae-ddf7e0f1b7ec" alt="C_all_tree7" width="839.7" height="750">

**Maximum-likelihood phylogenetic tree of the condensation domain superfamily.**

Subtype classification and sequences are described in the main text and the Method. Different subtypes are indicated by colors, with subtypes exclusive to fungi marked by underlines, and subtypes found predominantly in bacteria marked by asterisks. This tree is rooted, taking papA and WES as outgroups(black shading). L-clade and D-clade are indicated by blue and red shading, respectively.

**Note:** In the NRPS-motif-Finder result, E domain is not considered to be a kind of C domain subtype. And E domain has 7 motifs while C domain has 10 motifs.

# How to use
Please see [Wiki](https://github.com/RuolinHe/NRPS-motif-Finder/wiki)https://github.com/RuolinHe/NRPS-motif-Finder/wiki page.

# How to cite
If you have found antiSMASH useful, please [cite us](https://doi.org/10.1371/journal.pcbi.1011100).

# License
The analysis codes of this work are licensed under a [GNU General Public License-3.0 license](https://github.com/RuolinHe/NRPS-motif-Finder#GPL-3.0-1-ov-file).
