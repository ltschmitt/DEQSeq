#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

# read arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("OUTPREFIX, UMILOC must be provided", call.=FALSE)
}
OUTPREFIX = args[1]
UMILOC = args[2]

# split locations
UMILOC = unlist(strsplit(UMILOC,','))

# process each umi location (usually it should only be one)
umis = do.call(c,lapply(UMILOC, function(x){
      # extract umis from alignment
      astack = GenomicAlignments::stackStringsFromBam(paste0(OUTPREFIX,'1_preprocessing/regionfiltered.bam'), param=x, use.names = T)
      umis = as.character(astack) %>% gsub('-','',.,fixed = T)
      # filter umis by length
      umilen = as.integer(gsub('.+-','',x)) - as.integer(gsub('.+:|-[0-9]+$','',x)) + 1
      threshold = as.numeric(umilen)*0.8
      return(umis[str_length(umis) > threshold])
}))

# write
dir.create(paste0(OUTPREFIX,'2_3_vsearch'))
paste0('>',names(umis),'\n',umis) %>% writeLines(text = ., con = paste0(OUTPREFIX,'2_3_vsearch/umis.fasta'), sep = '\n')

