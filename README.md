# DEQSeq processing pipeline

# Requirements
Memory requirements: ~2x input fastq size
Hard drive requirements: ~7x input fastq size

Programs: R, samtools, minimap2, parallel, racon, medaka, vsearch
R packages: GenomicAlignments, tidyverse

# Installation

I recommend to install medaka over conda: 
```
conda create -n medaka -c conda-forge -c bioconda medaka
```
This has the advantage that you download and environment that already contains samtools, minimap2, medaka.
Then you only need to install racon, parallel, vsearch, and the R packages

```
conda install -c bioconda racon
git clone https://github.com/torognes/vsearch
```

Install R packages in R:
```{R}
# tidyverse installation
install.packages('tidyverse')
install.packages('stringdist')
install.packages('data.table')

# GenomicAlignments installation
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicAlignments")
```

Then download the pipeline
```
git clone https://github.com/ltschmitt/DEQSeq
```

# Usage
Fill out the yaml as described in the template and run the pipeline (in conda if you used it to install medaka):
```
conda activate medaka
bash DEQSeq/DEQSeq.sh -y template.yaml
```

The pipeline needs basecalled reads, everything with mean quality of Phred 10 or higher should suffice. You also need to prepare references for the Protein UMI region and the target site region. The references should be fairly similar to the sequenced reads, so that mapping reads to the reference with minimap2 is possible.

The pipeline expects no indels in the evolution products, insertions will be ignored and deletions will be filled with the reference sequence. This is necessary to ensure conversion to protein to work. If you expect indels, you can still use the pipeline without the last step, the polished reads contain all indels, but they are removed in the last step.

Edit calling can be done in two ways:
- Prepare references for each possible editing outcome, prefered method for large scale editing (e.g. recombination)
- Prepare one reference and define the location you want to extract for a detailed look, prefered method for base editing. Keep in mind that nanopore sequencing data can contain a lot of errors, so choosing a large region will yield a lot of different outcomes. Restricting the area to a few bp is therefore recommended.


