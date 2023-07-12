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

Parallel is often already installed on linux, if not please look up [GNU parallel](www.gnu.org/software/parallel) for how to install it.

```
# racon
conda install -c bioconda racon
# vsearch
git clone https://github.com/torognes/vsearch
cd vsearch
./autogen.sh
./configure CFLAGS="-03" CXXFLAGS="-03"
make
make install # as root or sudo make install
```

Install R packages in R:
```{R}
# tidyverse installation
install.packages('tidyverse')
install.packages('stringdist')

# GenomicAlignments installation
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicAlignments")
BiocManager::install("Biostrings")
```

Then download the pipeline
```
git clone https://github.com/ltschmitt/DEQSeq
```

# Usage

All basecalled sequence reads need to be in one sequence file (either .fastq or .fastq.gz; "cat" can be used for both of these files to combine the nanopore output into one file).

An optional preprocessing step can be performed by filtering for size and quality using [filtlong](https://github.com/rrwick/Filtlong).

You also need to prepare references for the Protein UMI region and the target site region (please see template.yaml for more information). The references should be fairly similar to the sequenced reads, so that mapping reads to the reference with minimap2 is possible.

Fill out the yaml as described in the template and run the pipeline (in conda if you used it to install medaka):
```
conda activate medaka
bash DEQSeq/DEQSeq.sh -y template.yaml
```

The pipeline expects no indels in the evolution products, insertions will be ignored and deletions will be filled with the reference sequence. This is necessary to ensure conversion to protein to work. If you expect indels, you can still use the pipeline, and retrieve the polished sequences with indels from the "polished_sequences.fastq" file. The final output file "Counts_Seqs.csv" contains the DNA sequence without insertions ("Sequence") and without insertions or deletions "Gapfix_sequence". The gaps are filled with the reference sequence bases and the number of gaps are noted in the "Gaps" column. The "ProteinSequence" contains the single letter amino acid sequence converted from "Gapfix_sequence".

Edit calling can be done in two ways:
- Prepare references for each possible editing outcome, prefered method for large scale editing (e.g. recombination)
- Prepare one reference and define the location you want to extract for a detailed look, prefered method for base editing. Keep in mind that nanopore sequencing data can contain a lot of errors, so choosing a large region will yield a lot of different outcomes. Restricting the area to a few bp is therefore recommended.
