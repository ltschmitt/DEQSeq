#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(GenomicAlignments))

# read arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("OUTPREFIX, PROTUMIREF, UMILOC and PROTLOC must be provided", call.=FALSE)
}
OUTPREFIX = args[1]
PROTUMIREF = args[2]
UMILOC = args[3]
PROTLOC = args[4]

########## Extract UMIs and Proteins ##############

# split locations
UMILOC = unlist(strsplit(UMILOC,','))
PROTLOC = unlist(strsplit(PROTLOC,','))

# get reference names
ref = Biostrings::readDNAStringSet(PROTUMIREF) %>% as.character()

# extract umis from polished alignment
umis = do.call(c,lapply(UMILOC, function(x){
      GenomicAlignments::stackStringsFromBam(paste0(OUTPREFIX,'polished_sequences.bam'), param=x, use.names = T) %>% as.character(.)
})) %>% gsub('-','',.,fixed = T) %>% enframe('Cluster','UMI')

# extract umis from polished alignment
genes = do.call(bind_rows,lapply(PROTLOC, function(x){
      seqs = GenomicAlignments::stackStringsFromBam(paste0(OUTPREFIX,'polished_sequences.bam'), param=x, use.names = T) %>% as.character(.)
      tibble(Cluster = names(seqs), Reference = gsub(':.*','',x), RefStart = as.integer(gsub('.+:|-[0-9]+','',x)), RefEnd = as.integer(gsub('.+-','',x)), Sequence = seqs)
}))

# fill gaps with reference sequence
filler = function(Sequence, RefSeq) {
      s = strsplit(Sequence,'')
      r = strsplit(RefSeq,'') %>% unlist()
      out = do.call(c,lapply(s, function(x){
	    ind = grep('-',x)
	    x[ind] = r[ind]
	    paste0(x,collapse='')
      }))
      return(out)
}
seqs_df = genes %>% mutate(Gaps = str_count(Sequence,'-'), RefSeq = str_sub(ref[Reference],RefStart,RefEnd), Gapfix_sequence = ifelse(Gaps>0, filler(Sequence,RefSeq), Sequence)) %>% select(-RefSeq)

# translate to protein
seqs_df = seqs_df %>% mutate(ProteinSequence = as.character( Biostrings::translate( Biostrings::DNAStringSet(Gapfix_sequence)))) %>% inner_join(umis,'Cluster')

########## Alignment Evaluation ##############

# Protein alignments evaluation
evals = read_tsv(paste0(OUTPREFIX,'5_cluster_polished_aligned/cluster_eval.tsv'), show_col_types = FALSE) %>% mutate(Protein_score = as.integer(gsub('AS:i:','',Protein_score,fixed = T)))
# plot alignment quality
pe = evals %>% group_by(Cluster) %>% mutate(Median_score = median(Protein_score)) %>% ungroup() %>% filter(Median_score %in% sort(unique(Median_score))[1:20]) %>% arrange(-Median_score) %>% mutate(Cluster = fct_inorder(as.factor(Cluster))) %>% ggplot(aes(x = Cluster, y = Protein_score)) + geom_boxplot() + labs(subtitle = 'Alignment scores of cluster sequences to polished Protein + UMI, higher is better') + coord_flip()
#save plot
nclusts = length(unique(evals$Cluster))
ggsave(paste0(OUTPREFIX,'worst_clusters.pdf'), pe, width = 8, height = 6, limitsize = F)

# TS alignment evaluation
tscounts = read_tsv(paste0(OUTPREFIX,'6_ts_aligned/ts_editing.tsv'), show_col_types = FALSE) %>% mutate(TS_score = as.integer(gsub('AS:i:','',TS_score,fixed = T)))
# plot alignment quality
pts = tscounts %>% ggplot(aes(x = TS, y = TS_score)) + geom_boxplot() + labs(subtitle = 'Alignment scores of TS alignments, higher is better')
#save plot
ggsave(paste0(OUTPREFIX,'ts_eval.pdf'), pts, width = 8, height = 6)

# combine and set alignment score QC outcome
tsce = inner_join(evals,tscounts, c('Cluster','ID')) %>% group_by(TS) %>% mutate(PassQC = TS_score > 0.75*median(TS_score) & Protein_score > 0.75*median(Protein_score))

# write
tsce %>% write_csv(paste0(OUTPREFIX,'Read_Stats.csv'))

########## Counting ##############

# filter for QC check and summarise/count
tsce_sum = tsce %>% filter(PassQC) %>% group_by(Cluster,TS) %>% summarise(n = n(), Mean_protein_score = round(mean(Protein_score),0), Mean_ts_score = round(mean(TS_score),0), .groups = 'drop')

# combine and write
inner_join(tsce_sum,seqs_df, 'Cluster') %>% write_csv(paste0(OUTPREFIX,'Counts_Seqs.csv'))

# Only when ts_location was specified for sequence extraction
# makes seperate count for all the detected ts sequences > 0.5 % occurance
if(suppressWarnings(length(tsce$TS_sequence)) > 0){
      tsce_seqsum = tsce %>% filter(PassQC) %>% group_by(Cluster,TS_sequence) %>% count() %>% group_by(Cluster) %>% mutate(Percent = n/sum(n)*100) %>% group_by(Cluster,TS_sequence)
      write_csv(tsce_seqsum,paste0(OUTPREFIX,'TS_sequence_counts.csv'))
}

