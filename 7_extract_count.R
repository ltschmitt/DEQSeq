#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(GenomicAlignments))

# read arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("OUTPREFIX, PROTUMIREF, UMILOC, PROTLOC, and MINCLUSTSIZE must be provided", call.=FALSE)
}
OUTPREFIX = args[1]
PROTUMIREF = args[2]
UMILOC = args[3]
PROTLOC = args[4]
MINCLUSTSIZE = args[5]

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

# combine and set alignment score QC outcome, I am penalizing bad protein scores more, a mix up of proteins is more likely than a messed up ts sequence
tsce = inner_join(evals,tscounts, c('Cluster','ID')) %>% group_by(TS) %>% mutate(PassQC = TS_score > 0.6*median(TS_score) & Protein_score > 0.8*median(Protein_score))
# TODO: multiple proteins might have multiple lengths - therefore I have to make a mean score for each protein when dealing with multiple proteins

# write
tsce %>% write_csv(paste0(OUTPREFIX,'Read_Stats.csv'))

########## Counting ##############

# filter reads and clusters for quality
mean_prot_score = tsce %>% pull(Protein_score) %>% mean()	
tsce = tsce %>% group_by(Cluster) %>% filter(mean(Protein_score) > mean_prot_score * 0.9) %>% filter(sum(PassQC)/n() > 0.8) %>% filter(PassQC)

# summarise/count
tsce_sum = tsce %>% group_by(Cluster,TS) %>% summarise(n = n(), Mean_protein_score = round(mean(Protein_score),0), Mean_ts_score = round(mean(TS_score),0), .groups = 'drop') %>% group_by(Cluster) #%>% filter(n >= MINCLUSTSIZE)

# combine
dat = inner_join(tsce_sum,seqs_df, 'Cluster') 

# join umis and filter by target read counts
UMI = dat %>% group_by(Cluster,UMI) %>% summarise(.groups = 'drop') %>% pull(UMI,Cluster)
udist = stringdist::stringdistmatrix(UMI, method = 'lv', useNames = 'names')
hc = hclust(udist)
umiclusters = cutree(hc, h = 3) %>% enframe(name = 'Cluster',value = 'NewCluster') 
udat = inner_join(umiclusters,dat, 'Cluster') %>% group_by(NewCluster,Reference,RefStart,RefEnd) %>% arrange(-n) %>% mutate(MainCluster = Cluster[1], SubClusters = paste0(unique(Cluster), collapse = ',')) %>% group_by(MainCluster,SubClusters,TS,Reference,RefStart,RefEnd) %>% summarise(n = sum(n), Gaps = Gaps[1], Sequence = Sequence[1], Gapfix_sequence = Gapfix_sequence[1], ProteinSequence = ProteinSequence[1], UMI = UMI[1],.groups = 'drop')
# subclusters are not the same if only some of them have reads for one specific target site, need to combine for all TS

# filter for min reads and write
udat %>% write_csv(paste0(OUTPREFIX,'Counts_Seqs.csv'))

# Only when ts_location was specified for sequence extraction
# makes seperate count for all the detected ts sequences > 0.5 % occurance
if(suppressWarnings(length(tsce$TS_sequence)) > 0){
      #tsce_seqsum = tsce %>% group_by(Cluster,TS_sequence) %>% count() %>% group_by(Cluster) %>% mutate(Percent = n/sum(n)*100) %>% group_by(Cluster,TS_sequence)
      tsce_seqsum = tsce %>% inner_join(umiclusters,'Cluster') %>% group_by(Cluster,TS) %>% arrange(-n()) %>% group_by(NewCluster,TS) %>% mutate(MainCluster = Cluster[1]) %>% group_by(MainCluster,TS_sequence) %>% count() %>% group_by(MainCluster) %>% mutate(Percent = n/sum(n)*100) %>% group_by(MainCluster,TS_sequence)
      write_csv(tsce_seqsum,paste0(OUTPREFIX,'TS_sequence_counts.csv'))
}

# TODO: check location ts umi joined
