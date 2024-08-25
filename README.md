# zAMExplorer
zAMPExplorer is an interactive Shiny application designed to facilitate the downstream analysis of outputs generated by the [zAMP pipeline](https://zamp.readthedocs.io/en/latest/) . zAMP itself is a comprehensive DADA2-based bioinformatics pipeline for processing 16S rRNA gene amplicon metagenomic sequences, offering convenient, reproducible, and scalable analysis from raw sequencing reads to taxonomic profiles. The output of zAMP is a phyloseq object, which serves as the input for zAMPExplorer.

zAMPExplorer enables users to perform a wide range of microbiota and statistical analyses, including compositional barplot, relative abundance heatmap, community diversity (alpha diversity), community similarity through unsupervised (NMDS/PCoA) and supervised (RDA) ordinations, and community typing (or clustering) of microbial profiles using Dirichlet Multinomial Mixtures (DMM). All of these analyses are made accessible through an intuitive graphical interface, bridging the gap between complex command-line bioinformatics processing and user-friendly data exploration.
