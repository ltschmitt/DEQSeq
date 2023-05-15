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
ucf = uc %>% group_by(Cluster) %>% mutate(Size = n()) %>% filter(Size >= as.integer(MINCLUSTSIZE)) %>% arrange(-Size,Cluster) %>% mutate(Cluster = as.integer(fct_inorder(as.factor(as.character(Cluster))))) 

write(paste0('Clusters: ', length(unique(ucf$Cluster)), '\n', 'Reads used: ', length((ucf$Cluster))), file = paste0(OUTPREFIX,'cluster_stats.txt'))
message(paste0('Clusters: ', length(unique(ucf$Cluster)), '    ', 'Reads used: ', length((ucf$Cluster))))

# make folder to safe split reads
dir.create(paste0(OUTPREFIX,'4_cluster_fastq'))

# open file for chunk reading
con = file(paste0(OUTPREFIX,"1_preprocessing/regionfiltered.tsv"), 'r') # TODO: change this to read the fastq so no tsv has to be made
chunksize = 10000 # determines how many lines are read at once
# loop through the lines of the file, combine with cluster table (ucf) and write into seperate files
while(T) {
      chunk = readLines(con, chunksize)
      if(length(chunk) == 0) break
      else {
	    chunk %>% as_tibble() %>% separate_wider_delim(value, delim = '\t', names = c("ID","Seq", "Qual")) %>% inner_join(ucf,'ID') %>% transmute(Cluster, output = paste0('@',ID,'\n',Seq,'\n+\n',Qual)) %>% group_by(Cluster) %>% group_walk(~ write(.x$output,file = paste0(OUTPREFIX,'4_cluster_fastq_test/Cluster_',unique(.y$Cluster),'.fastq'), sep = '\n', append = T))
      }
}

# TODO: check if files are already present before append writing the split fastqs.

# data table variant, that uses more memory
#data.table::fread(cmd = paste0("samtools view ",OUTPREFIX,"1_preprocessing/regionfiltered.bam | awk '{print $1,$10,$11}'"), col.names = c('ID','Seq','Qual'), quote = "") %>% inner_join(ucf,'ID') %>% transmute(Cluster, output = paste0('@',ID,'\n',Seq,'\n+\n',Qual)) %>% group_by(Cluster) %>% group_walk(~ write(.x$output,file = paste0(OUTPREFIX,'4_cluster_fastq_test/Cluster_',unique(.y$Cluster),'.fastq'), sep = '\n'))

