crAssphage project
================

-   [Background](#background)
    -   [Introduction](#introduction)
-   [Bioinformatics](#bioinformatics)
    -   [Phage genomes and Bowtie index](#phage-genomes-and-bowtie-index)
    -   [Mapping reads against phage genomes and calculating genome coverage](#mapping-reads-against-phage-genomes-and-calculating-genome-coverage)
    -   [ARG and MGE abundance](#arg-and-mge-abundance)
-   [Data analysis and statistics in R](#data-analysis-and-statistics-in-r)
    -   [Load in the data and libraries needed in the analyses](#load-in-the-data-and-libraries-needed-in-the-analyses)
    -   [Figure 1 - crAssphage in human fecal metagenomes](#figure-1---crassphage-in-human-fecal-metagenomes)
    -   [Figure 2 - crAssphage in impacted environments](#figure-2---crassphage-in-impacted-environments)
    -   [Figure 3 - crAsspahge in MG-RAST metagenomes](#figure-3---crasspahge-in-mg-rast-metagenomes)
-   [References](#references)

Background
==========

Supplementary data analysis for the paper:

***Karkman, A., Pärnänen, K., Larsson, DGJ.*** 2018. Fecal pollution explains antibiotic resistance gene abundances in anthropogenically impacted environments. *Nature Communications* X: XYZ. <https://doi.org/XYZ>.

Introduction
------------

Discharge of treated sewage introduces antibiotic resistance genes (ARGs) and resistant bacteria (ARBs) to the environment (Karkman et al. 2018). It has been speculated that the resistance determinants could be selected and/or disseminated to other bacteria in the downstream environments due to the antibiotcs and other selective agents released in the effluents (Rizzo et al. 2013; Guo et al. 2017). However, the increased abundance of ARGs in downstream environments could as well be explained by the amount of fecal pollution without any large scale selection and/or dissemination of the genes.

Recently a human feacal phage, crAssphage, was discovered from human fecal metagenomes (Dutilh et al. 2014) and was shown to infect *Bacteroides intestinalis* (Shkoporov et al. 2018). It has been shown to be very abundant in human fecal material and mostly specific to humans (García-Aljaro et al. 2017). Due to these facts, it has already been used as fecal marker in various studies (Stachler et al. 2018; Ahmed, Zhang, et al. 2018; Stachler and Bibby 2014; Stachler et al. 2017; Ahmed, Lobos, et al. 2018). Another fecal Bacteroides phage, ɸB124-14, has also been used for microbial source tracking. Unlike crAssphage, it should be abundant in porcine and bovine guts as well (Ogilvie et al. 2017). However, we did not find it to perform as well as crAssphage and it won't be included in the analyses.

In this study we show that in most of the the studied environments the abundance of ARGs, *intI1* integrase gene and mobile genetic elements correlates well with fecal pollution levels with no evident signs of selection or dissemination of the resistance genes. The only exception being sediments polluted with wastewater from drug manufacturing containing exceptionally high levels of antibiotics (Bengtsson-Palme et al. 2014; Kristiansson et al. 2011).

Bioinformatics
==============

The bioinformartics part shows only example commands and is not meant to be run as such. The results from the bioinformatics paret are available in the `data` folder and will be used in the datraanalysis part with R.

All metagenomic samples were downloaded from public repositories as described in the methods.

The crAssphage ([NC\_024711.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_024711.1)) and ɸB124-14 ([HE608841.1](https://www.ncbi.nlm.nih.gov/nuccore/HE608841.1)) genomes were downloaded from GenBank as fasta files.

Phage genomes and Bowtie index
------------------------------

The phage genomes were indexed using `bowtie2-build`

``` bash
bowtie2-build phage_genome.fasta phage_genome
```

Mapping reads against phage genomes and calculating genome coverage
-------------------------------------------------------------------

After indexing the phage genomes each sample was mapped against the genomes. The average genome coverage was calculated and used as a proxy for phage abundance.
**Paired-end reads:**

``` bash
bowtie2 -x phage_genome -1 Sample_R1_reads.fastq.gz -2 Sample_R2_reads.fastq.gz -S Sample.sam
samtools view -Sb -f 2 Sample.sam > Sample.bam
samtools sort Sample.bam -o Sample_sort.bam
samtools index Sample_sort.bam
export GEN_COV=$(samtools depth -a Sample_sort.bam |\
                  awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
echo 'Sample\t'$GEN_COV
```

**Single reads:**

``` bash
bowtie2 -x phage_genome -U Sample.fastq -S Sample.sam
samtools view -Sb -q 10 Sample.sam > Sample.bam
samtools sort  Sample.bam -o  Sample_sort.bam
samtools index Sample_sort.bam
export GEN_COV=$(samtools depth -a Sample_sort.bam |\
                  awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
echo 'Sample\t'$GEN_COV 
```

ARG and MGE abundance
---------------------

The `fastq` files were converted to `fasta` and the sample name was added to each sequence header before concatenating all fasta files to one file for annotation with DIAMOND.

The count tables are generated from the DIAMOND outputs using custom scripts. The scripts for single (`parse_diamond.py`) and paired-end reads (`parse_diamondPE.py`) can be downloaded from [here](https://github.com/karkman/parse_diamond).

**Paired-end reads:**

``` bash
diamond blastx -d ResFinder -q Sample_R1.fasta --max-target-seqs 1 -o Sample_R1_res.txt \
                -f 6 --id 90 --min-orf 20 -p 24 --masking 0
diamond blastx -d ResFinder -q Sample_R2.fasta --max-target-seqs 1 -o Sample_R2_res.txt \
                -f 6 --id 90 --min-orf 20 -p 24 --masking 0
python parse_diamondPE.py -1 Sample_R1_res.txt -2 Sample_R2_res.txt -o Sample_ResFinder.csv
```

**Single reads:**

``` bash
diamond blastx -d ResFinder -q Sample.fasta --max-target-seqs 1 -o Sample_res.txt \
                -f 6 --id 90 --min-orf 20 -p 24 --masking 0
python parse_diamondPE.py -i Sample_res.txt -o Sample_ResFinder.csv
```

To bp count in each metagenome was used for normalization.

Data analysis and statistics in R
=================================

The results from mapping against crAssphage and gene annotations were imported to R.

Load in the data and libraries needed in the analyses
-----------------------------------------------------

All the data can be found from the `data` folder.

``` r
library(tidyverse)
HMP <- read.table("data/")
crass_imapct <- read.table("data/")
crass_wwtp <- read.table("data/")
crass_MG_RAST <- read.table("data/")
res_risk <- read.table("data/")
```

Figure 1 - crAssphage in human fecal metagenomes
------------------------------------------------

text here

``` r
par(fig=c(0,0.45,0,0.8), new=TRUE)
plot(log10(rel_res)~log10(rel_crAss), data=HMP, bg=cols[as.factor(HMP$country)], pch=21,
      ylab = "Normalized ARG abundance (log10)", xlab="Normalized crAssphage abundance (log10)", cex=2, ylim=c(2.5, 4.5))
par(fig=c(0,0.45,0.5,1), new=TRUE)
boxplot(log10(rel_crAss)~country, data=HMP, horizontal=TRUE, col=cols, axes=F)
axis(2, at=1:3, labels=c("China", "Europe", "US"), las=1)
title("A", adj = 0, line = 0)

par(fig=c(0.45,0.9,0,0.8), new=TRUE)
tmp <- subset(HMP, rel_int>0)
plot(log10(rel_res)~log10(rel_int), data=tmp, bg=cols[as.factor(tmp$country)], pch=21,
      ylab = "", xlab="Normalized intI1 abundance (log10)", cex=2, ylim=c(2.5, 4.5))
par(fig=c(0.45,0.9,0.5,1), new=TRUE)
boxplot(log10(rel_int)~country, data=tmp, horizontal=TRUE, col=cols, axes=F)
axis(2, at=1:3, labels=c("China", "Europe", "US"), las=1)
title("B", adj = 0, line = 0)

par(fig=c(0.8,1,0,0.8),new=TRUE)
boxplot(log10(rel_res)~country, data=HMP, col=cols, axes=F)
axis(1, at=1:3, labels=c("China", "Europe", "US"), las=3)
```

Figure 2 - crAssphage in impacted environments
----------------------------------------------

text here

``` r
ggplot(crass_impact, aes(x=rel_crAss, y=rel_res, color=country)) +  geom_smooth(method="lm") + geom_point(aes(shape=crAss_detection), size=5) + scale_x_log10() + scale_y_log10() + theme_classic() +
            labs(y = "Normalized ARG abundance", x="Normalized crAssphage abundance", color="Study", shape="crAssphage detection") + scale_colour_manual(values=cols)
```

Figure 3 - crAsspahge in MG-RAST metagenomes
--------------------------------------------

text here.

``` r
ggplot(MG_RAST_crass, aes(x=rel_crAss, y=rel_res, color=revised_source)) + geom_point(size=5) + scale_x_log10() + scale_y_log10() + geom_smooth(method="lm") + theme_classic() +
      labs(y = "Normalized ARG abundance", x="Normalized crAssphage abundance", color = "Revised source") + scale_colour_manual(values=cols)
```

References
==========

Ahmed, Warish, Aldo Lobos, Jacob Senkbeil, Jayme Peraud, Javier Gallard, and Valerie J. Harwood. 2018. “Evaluation of the Novel crAssphage Marker for Sewage Pollution Tracking in Storm Drain Outfalls in Tampa, Florida.” *Water Research* 131 (March): 142–50. <https://doi.org/10.1016/j.watres.2017.12.011>.

Ahmed, Warish, Qian Zhang, Aldo Lobos, Jacob Senkbeil, Michael J. Sadowsky, Valerie J. Harwood, Nazanin Saeidi, Oswald Marinoni, and Satoshi Ishii. 2018. “Precipitation Influences Pathogenic Bacteria and Antibiotic Resistance Gene Abundance in Storm Drain Outfalls in Coastal Sub-Tropical Waters.” *Environment International* 116 (July): 308–18. <https://doi.org/10.1016/j.envint.2018.04.005>.

Bengtsson-Palme, Johan, Fredrik Boulund, Jerker Fick, Erik Kristiansson, and D. G. Joakim Larsson. 2014. “Shotgun Metagenomics Reveals a Wide Array of Antibiotic Resistance Genes and Mobile Elements in a Polluted Lake in India.” *Frontiers in Microbiology* 5 (December): 648. <https://doi.org/10.3389/fmicb.2014.00648>.

Dutilh, Bas E., Noriko Cassman, Katelyn McNair, Savannah E. Sanchez, Genivaldo G. Z. Silva, Lance Boling, Jeremy J. Barr, et al. 2014. “A Highly Abundant Bacteriophage Discovered in the Unknown Sequences of Human Faecal Metagenomes.” *Nature Communications* 5: 4498. <https://doi.org/10.1038/ncomms5498>.

García-Aljaro, Cristina, Elisenda Ballesté, Maite Muniesa, and Juan Jofre. 2017. “Determination of crAssphage in Water Samples and Applicability for Tracking Human Faecal Pollution.” *Microbial Biotechnology* 10 (6): 1775–80. <https://doi.org/10.1111/1751-7915.12841>.

Guo, Jianhua, Jie Li, Hui Chen, Philip L. Bond, and Zhiguo Yuan. 2017. “Metagenomic Analysis Reveals Wastewater Treatment Plants as Hotspots of Antibiotic Resistance Genes and Mobile Genetic Elements.” *Water Research* 123 (October): 468–78. <https://doi.org/10.1016/j.watres.2017.07.002>.

Karkman, Antti, Thi Thuy Do, Fiona Walsh, and Marko P.J. Virta. 2018. “Antibiotic-Resistance Genes in Waste Water.” *Trends in Microbiology* 26 (3): 220–28. <https://doi.org/10.1016/j.tim.2017.09.005>.

Kristiansson, Erik, Jerker Fick, Anders Janzon, Roman Grabic, Carolin Rutgersson, Birgitta Weijdegård, Hanna Söderström, and D. G. Joakim Larsson. 2011. “Pyrosequencing of Antibiotic-Contaminated River Sediments Reveals High Levels of Resistance and Gene Transfer Elements.” Edited by Francisco Rodriguez-Valera. *PLoS ONE* 6 (2): e17038. <https://doi.org/10.1371/journal.pone.0017038>.

Ogilvie, Lesley A., Jonathan Nzakizwanayo, Fergus M. Guppy, Cinzia Dedi, David Diston, Huw Taylor, James Ebdon, and Brian V. Jones. 2017. “Resolution of Habitat-Associated Ecogenomic Signatures in Bacteriophage Genomes and Application to Microbial Source Tracking.” *The ISME Journal*, December. <https://doi.org/10.1038/s41396-017-0015-7>.

Rizzo, L., C. Manaia, C. Merlin, T. Schwartz, C. Dagot, M.C. Ploy, I. Michael, and D. Fatta-Kassinos. 2013. “Urban Wastewater Treatment Plants as Hotspots for Antibiotic Resistant Bacteria and Genes Spread into the Environment: A Review.” *Science of the Total Environment* 447 (March): 345–60. <https://doi.org/10.1016/j.scitotenv.2013.01.032>.

Shkoporov, Andrey N., Ekaterina V. Khokhlova, C. Brian Fitzgerald, Stephen R. Stockdale, Lorraine A. Draper, R. Paul Ross, and Colin Hill. 2018. “ΦCrAss001 Represents the Most Abundant Bacteriophage Family in the Human Gut and Infects Bacteroides Intestinalis.” *Nature Communications* 9 (1): 4781. <https://doi.org/10.1038/s41467-018-07225-7>.

Stachler, Elyse, Benay Akyon, Nathalia Aquino de Carvalho, Christian Ference, and Kyle Bibby. 2018. “Correlation of crAssphage qPCR Markers with Culturable and Molecular Indicators of Human Fecal Pollution in an Impacted Urban Watershed.” *Environmental Science & Technology* 52 (13): 7505–12. <https://doi.org/10.1021/acs.est.8b00638>.

Stachler, Elyse, and Kyle Bibby. 2014. “Metagenomic Evaluation of the Highly Abundant Human Gut Bacteriophage CrAssphage for Source Tracking of Human Fecal Pollution.” *Environmental Science & Technology Letters* 1 (10): 405–9. <https://doi.org/10.1021/ez500266s>.

Stachler, Elyse, Catherine Kelty, Mano Sivaganesan, Xiang Li, Kyle Bibby, and Orin C. Shanks. 2017. “Quantitative CrAssphage PCR Assays for Human Fecal Pollution Measurement.” *Environmental Science & Technology* 51 (16): 9146–54. <https://doi.org/10.1021/acs.est.7b02703>.
