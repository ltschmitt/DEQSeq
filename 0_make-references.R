#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

# read arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("OUTPREFIX, PROTUMIREF must be provided", call.=FALSE)
}
OUTPREFIX = args[1]
PROTUMIREF = args[2]
TSREF = args[3]

# Function to read fasta files
read_fasta = function(input) {
      dat = readLines(input)
      dat = dat[dat != '']
	    seqs = grep('^[^>]',dat)
	    seqs = tibble::tibble(Group = cumsum(c(1,diff(seqs)) !=1), Sequence = dat[seqs]) %>% dplyr::group_by(Group) %>% dplyr::summarise(Sequence = paste0(Sequence,collapse = '')) %>% dplyr::pull(Sequence)
	    return(stats::setNames(seqs, gsub('>', '', grep('^>',dat,value=T), fixed = T)))
}

dir.create(paste0(OUTPREFIX,'0_references'), recursive = T)

# PROTUMIREF bed files
protref = read_fasta(PROTUMIREF)
tibble(Name = names(protref), C1 = 30, C2 = 31) %>% write_tsv(paste0(OUTPREFIX,'0_references/protumi1.bed'), col_names = FALSE)
tibble(Name = names(protref), C1 = str_length(protref)-30, C2 = str_length(protref)-29) %>% write_tsv(paste0(OUTPREFIX,'0_references/protumi2.bed'), col_names = FALSE)

# TSREF bed files
tsref = read_fasta(TSREF)
tibble(Name = names(tsref), C1 = 30, C2 = 31) %>% write_tsv(paste0(OUTPREFIX,'0_references/ts1.bed'), col_names = FALSE)
tibble(Name = names(tsref), C1 = str_length(tsref)-30, C2 = str_length(tsref)-29) %>% write_tsv(paste0(OUTPREFIX,'0_references/ts2.bed'), col_names = FALSE)
