<img src="https://github.com/RuolinHe/NRPS-motif-Finder/assets/76482251/5ee2d1fc-b026-4be3-bfdc-6ac6d6151a8b" alt="Striking Image" width="600" height="600">

The [cover](https://journals.plos.org/ploscompbiol/issue?id=10.1371/issue.pcbi.v19.i05) of PLOS Computational Biology in 2023 May for this work.

# What is NRPS-motif-Finder
NRPS-motif-Finder is a tool for standardization of Non-ribosomal peptide synthetase (NRPS). It partitions the input NRPS protein sequence by locating these conserved motif, to output a motif-and-intermotif architecture that feeds in subsequent analysis such as C domain classification, NRPS re-engineering.

<img src="https://github.com/RuolinHe/NRPS-motif-Finder-matlab-version/assets/76482251/a0358bd6-3431-4bda-9477-20898debc41a" alt="logo" width="388.9" height="200">

***

**Note**: NRPS-motif-Finder is **the only one** that can classify **all fungal** C domain subtypes so far!

***

# Supported domains and motifs
Adenylation (A) domain has 12 domain: Aalpha, A1-A5, G-motif, A6-A10. Among them, Aalpha and G-motif were two new motifs proposed in [our paper](https://doi.org/10.1371/journal.pcbi.1011100).

Condensation (C) domain has 10 domain: C1-C10.

Thiolation (T) domain has 2 domain: Talpha, T1. Talpha was one new motif [our paper](https://doi.org/10.1371/journal.pcbi.1011100).

Thioesterase (TE) domain has 1 domain: TE1.

Epimerization (E) domain has 7 domains: E1-E7.

# Featured function
For more details about featured functions, please check [our Wiki](https://github.com/RuolinHe/NRPS-motif-Finder/wiki/1%E2%80%90Featured%E2%80%90function).
## C domain subtype classification
One of the most important features is that NRPS-motif-Finder supports the full subtype classification of C domain.

<img src="https://github.com/RuolinHe/NRPS-motif-Finder/assets/76482251/6fef2845-67aa-44f6-8f7d-7638805b45c9" alt="C_all_tree7" width="839.7" height="750">

**Maximum-likelihood phylogenetic tree of the condensation domain superfamily.**

Subtype classification and sequences are described in the main text and the Method. Different subtypes are indicated by colors, with subtypes exclusive to fungi marked by underlines, and subtypes found predominantly in bacteria marked by asterisks. This tree is rooted, taking papA and WES as outgroups(black shading). L-clade and D-clade are indicated by blue and red shading, respectively.

**Note:** In the NRPS-motif-Finder result, E domain is not considered to be a kind of C domain subtype. And E domain has 7 motifs while C domain has 10 motifs.

## A domain loop group classification
Loop length and loop group are proposed in [our paper](https://doi.org/10.1371/journal.pcbi.1011100). And we found loop group is related with A domain substrate specificity.

![Fig3_20230113](https://github.com/RuolinHe/NRPS-motif-Finder/assets/76482251/b340f9fd-a0ae-42e9-9c8c-ecda58ac3c18)

**The specificity-conferring code of the A domain is correlated with loop length and phylogeny**

**A.** SCA of 2,636 A domain sequences, together with their substrate specificities attached to the last column of the multiple sequence alignment. Six sectors with a high contribution from the substrate column (>0.05, the size of points on the left scales the substrate’s contribution to the sector, see Method for details) are sorted by their eigenvalues. The size of points scales its contribution to the sector. Orange bars mark the A domain motifs from A1 to A8. The start and end of the five loop regions are marked by black and green dotted lines, respectively. S4 and S6 are the 4th and 6th of the specificity-conferring codes. G is the G-motif.

**B.** Distance matrix of A domain. Upper right on the heatmap is the Euclidean distance of the loop length as a 5-element vector. Lower left on the heatmap is the sequence distance of the A domain. The matrix is sorted by the substrate specificity followed by the loop length group. Substrates, groups of loop length, and phylum of these A domains, are shown by colors in sidebars.

**C.** Example showing that A domains conferring identical substrate exhibit distinct specificity-conferring codes, when they are categorized into different loop-length-groups. Phylum composition in each group is shown in the pie chart.

## New motifs
We found there are some conserved sites in A domain and T domain, and they weren't known motif.

The new motif before A1 was named as Aalpha motif. The new motif between A5 and A6 was named as G-motif due to the core conserved Gly in the motif. The new motif before T1 was named as Talpha motif.

<img src="https://github.com/RuolinHe/NRPS-motif-Finder/assets/76482251/89244090-e159-4eab-af7c-2561e83a823c" alt="Fig5_revision20230113A" width="800" height="586.8">

**Analysis of amino acid frequency reveals potential new motifs**

Amino acid frequency and gap frequency along the multiple sequence alignment of the NRPS CAT modules. In the bottom panel, bar heights indicate the frequency of the most frequent amino acid. Bars in the known core motifs from the C, A, and T domains were colored blue, orange, and yellow, respectively. The horizontal red dashed line represents the 0.95 frequency level. Domain boundaries annotated by Pfam are divided by red triangles. The colored patch above the amino acid frequency indicates gap frequency. Three potential new motifs (position 1183, 1960, and 2435/2447 in MSA) are marked by the blue dashed box. The upper panel shows the sequence logo and the gap frequency near the three potential new motifs.

## NRPS motif sequence logo repository
In [our paper](https://doi.org/10.1371/journal.pcbi.1011100), we construct the most comprehensive NRPS motif sequence logo in 16,820 bacterial and 2,505 fungal genomes (as of 2022/10/23).

In this large dataset, we analyzed 83,489 C domains, 95,582 A domains, 86,688 T domains, 14,502 E domains, and 23,590 TE domains from bacteria; and 34,269 C domains, 40,458 A domains, 26,651 T domains, 3,982 E domains, and 4,008 TE domains from fungi.

<img src="https://github.com/RuolinHe/NRPS-motif-Finder/assets/76482251/42dd953a-a5a4-4950-adbb-6bc0836128cc" alt="C_subtype_figure_D" width="707.76" height="800">

**Figure 2D.** The sequence logo for the C3 or E2 motif from different C domain subtypes and the T1 or ACP1 motif adjacent to each subtype. Sequences from bacteria were marked by red, while sequences from fungi were marked by blue.

# How to use
Please see [Wiki](https://github.com/RuolinHe/NRPS-motif-Finder/wiki) page.

# How to cite
If you have found antiSMASH useful, please [cite us](https://doi.org/10.1371/journal.pcbi.1011100).

# License
The analysis codes of this work are licensed under a [GNU General Public License-3.0 license](https://github.com/RuolinHe/NRPS-motif-Finder#GPL-3.0-1-ov-file).
