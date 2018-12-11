crAssphage project
================

-   [Background](#background)
-   [Bioinformatics](#bioinformatics)
    -   [Phage genomes and Bowtie index](#phage-genomes-and-bowtie-index)
    -   [Mapping reads against phage genomes and calculating genome coverage](#mapping-reads-against-phage-genomes-and-calculating-genome-coverage)
    -   [ARG and MGE abundance](#arg-and-mge-abundance)
-   [Data analysis and statistics in R](#data-analysis-and-statistics-in-r)
    -   [Load in the data and libraries needed in the analyses](#load-in-the-data-and-libraries-needed-in-the-analyses)
    -   [Figure 1 - crAssphage in human fecal metagenomes](#figure-1---crassphage-in-human-fecal-metagenomes)

Background
==========

*Fecal pollution explains antibiotic resistance gene abundances in anthropogenically impacted environments*

The source code for all the analyses in the above paper can be found from here.

Bioinformatics
==============

Some text here

Phage genomes and Bowtie index
------------------------------

Text here

``` bash
bowtie2-build phage_genome.fasta phage_genome
```

Mapping reads against phage genomes and calculating genome coverage
-------------------------------------------------------------------

Mapping against phage genomes and calculation of the genome coverage as a proxy for phage abundance

``` bash
bowtie2 -x phage_genome  -p 2 -1 Sample_R1_reads.fastq.gz -2 Sample_R2_reads.fastq.gz -S Sample.sam
samtools view -Sb -f 2 Sample.sam > Sample.bam
samtools sort Sample.bam -o Sample_sort.bam
samtools index Sample_sort.bam
export GEN_COV=$(samtools depth -a Sample_sort.bam |awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
echo 'Sample\t'$GEN_COV
```

ARG and MGE abundance
---------------------

some text here as well. Link to the `parse_diamondPE.py` script. Link to Resfinder and MGE database.

``` bash
diamond blastx -d ResFinder -q Sample_R1.fasta --max-target-seqs 1 -o Sample_R1_res.txt \
                -f 6 --id 90 --min-orf 20 -p 24 --masking 0
diamond blastx -d ResFinder -q Sample_R2.fasta --max-target-seqs 1 -o Sample_R2_res.txt \
                -f 6 --id 90 --min-orf 20 -p 24 --masking 0
python parse_diamondPE.py -1 Sample_R1_res.txt -2 Sample_R2_res.txt -o Sample_ResFinder.csv
```

Data analysis and statistics in R
=================================

Some text here.

Load in the data and libraries needed in the analyses
-----------------------------------------------------

``` r
library(tidyverse)
HMP <- read.table()
crass_imapct <- read.table()
crass_wwtp <- read.table()
crass_MG_RAST <- read.table()
res_risk <- read.table()
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
