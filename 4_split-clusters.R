suppressPackageStartupMessages(library(tidyverse))

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("OUTPREFIX, MINCLUSTSIZE and UMIIDENT2 must be provided", call.=FALSE)
}
OUTPREFIX = args[1]
MINCLUSTSIZE = args[2]
UMIIDENT2 = suppressWarnings(as.numeric(args[3]))

# read clusters
hnames = c('Type','Cluster','length_size','Similarity','Orientation','x1','x2','CIGAR','ID','Centroid_ID')
uc = read_tsv(paste0(OUTPREFIX,'2_3_vsearch/uc.tab'), col_names = hnames,show_col_types = FALSE)[c(9,1:3)]
uc = uc %>% filter(Type %in% c('S','H'))

if(!is.na(UMIIDENT2)){
      uc2 = read_tsv(paste0(OUTPREFIX,'2_3_vsearch/uc2.tab'), col_names = hnames,show_col_types = FALSE)[c(9,1:3)]
      uc2 = uc2 %>% filter(Type %in% c('S','H')) %>% mutate(NewCluster = Cluster, Cluster = as.integer(gsub('.*clusterid=','',ID))) %>% select(Cluster,NewCluster)
      uc = inner_join(uc,uc2,'Cluster') %>% mutate(Cluster = NewCluster) %>% select(-NewCluster)
} 

# filter clusters
ucf = uc %>% group_by(Cluster) %>% mutate(Size = n()) %>% filter(Size >= as.integer(MINCLUSTSIZE))
write(paste0('Clusters: ', length(unique(ucf$Cluster)), '\n', 'Reads used: ', length((ucf$Cluster))), file = paste0(OUTPREFIX,'cluster_stats.txt'))
message(paste0('Clusters: ', length(unique(ucf$Cluster)), '    ', 'Reads used: ', length((ucf$Cluster))))

# read bam, renumber clusters and prep for writing
bam = data.table::fread(cmd = paste0("samtools view ",OUTPREFIX,"1_preprocessing/regionfiltered.bam | awk '{print $1,$10,$11}'"), col.names = c('ID','Seq','Qual'), quote = "")
out = bam %>% inner_join(ucf,'ID') %>% arrange(-Size,Cluster) %>% mutate(Cluster = as.integer(fct_inorder(as.factor(as.character(Cluster))))) %>% transmute(Cluster, output = paste0('@',ID,'\n',Seq,'\n+\n',Qual))

# write seperate files
dir.create(paste0(OUTPREFIX,'4_cluster_fastq'))
invisible(lapply(split(out,out$Cluster), function(x) write(x$output,file = paste0(OUTPREFIX,'4_cluster_fastq/Cluster_',unique(x$Cluster),'.fastq'), sep = '\n')))
