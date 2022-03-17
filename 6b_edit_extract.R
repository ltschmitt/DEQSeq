#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(tidyverse))

# read arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("OUTPREFIX, TSREF and TSLOC must be provided", call.=FALSE)
}
OUTPREFIX = args[1]
TSREF = args[2]
TSLOC = args[3]

########## Extract UMIs and Proteins ##############

# split locations
TSLOC = unlist(strsplit(TSLOC,','))

# get cluster names
clusters = list.files(paste0(OUTPREFIX,'6_ts_aligned/'), pattern = '.bam$')

# extract ts sequences from ts alignment

tsseqs = do.call(bind_rows,lapply(TSLOC, function(x){
      do.call(bind_rows,lapply(clusters, function(y){
	    astack = GenomicAlignments::stackStringsFromBam(paste0(OUTPREFIX,'6_ts_aligned/',y), param=x, use.names = T) %>% as.character(.)
	    return(tibble(ID = names(astack),TS_location = x, TS_sequence = astack))
      }))
})) %>% mutate(TS = gsub(':.+','',TS_location), TS_location = gsub('^.+:','',TS_location))

# read ts_editing file for Alignment score
dat = read_tsv(paste0(OUTPREFIX,'6_ts_aligned/ts_editing.tsv'), show_col_types = F) %>% inner_join(tsseqs,c('ID','TS'))

# overwrite ts_editing file, adds extracted TS site
write_tsv(dat, paste0(OUTPREFIX,'6_ts_aligned/ts_editing.tsv'))
