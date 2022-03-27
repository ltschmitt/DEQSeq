#!/bin/bash

#source /opt/miniconda3/etc/profile.d/conda.sh
#conda activate medaka141

OUTPREFIX=$1
PROTUMIREF=$2 
THREADS=$3
MODEL=$4
minimap2LOC=$5
samtoolsLOC=$6
medaka_consensusLOC=$7
raconLOC=$8
parallelLOC=$9

mkdir ${OUTPREFIX}5_cluster_aligned
mkdir ${OUTPREFIX}5_racon
mkdir ${OUTPREFIX}5_medaka
mkdir ${OUTPREFIX}5_cluster_polished_aligned

# prepare function for parallel processing
dowork() {
      file=$1
      PROTUMIREF=$2
      OUTPREFIX=$3
      MODEL=$4
      minimap2LOC=$5
      samtoolsLOC=$6
      medaka_consensusLOC=$7
      raconLOC=$8

      bname=$(basename -s .fastq $file)

      # align clusters to reference
      ${minimap2LOC} -ax map-ont -t 1 --secondary=no $PROTUMIREF $file > ${OUTPREFIX}5_cluster_aligned/$bname.unfilt.sam
      # reduce reads to most mapped reference if there are multiple references
      if [[ $(grep '^>' $PROTUMIREF | wc -l) -gt 1 ]] ; then
	    refname=$(${samtoolsLOC} view ${OUTPREFIX}5_cluster_aligned/$bname.unfilt.sam | cut -f 3 | sort | uniq -c | sort -n | tail -n 1 | awk '{print $2}')
	    ${samtoolsLOC} view -H ${OUTPREFIX}5_cluster_aligned/$bname.unfilt.sam > ${OUTPREFIX}5_cluster_aligned/$bname.sam
	    ${samtoolsLOC} view ${OUTPREFIX}5_cluster_aligned/$bname.unfilt.sam | awk -v refname=$refname '{ if ($3 == refname) { print } }' >> ${OUTPREFIX}5_cluster_aligned/$bname.sam
	    ${samtoolsLOC} fastq ${OUTPREFIX}5_cluster_aligned/$bname.sam > $file
	    #rm ${OUTPREFIX}5_cluster_aligned/$bname.unfilt.sam
      else
	    mv ${OUTPREFIX}5_cluster_aligned/$bname.unfilt.sam ${OUTPREFIX}5_cluster_aligned/$bname.sam
      fi
      # index fastq for racon and medaka
      ${samtoolsLOC} faidx $file
      # polish with racon
      ${raconLOC} -m 8 -x -6 -g -8 -w 500 --no-trimming -t 1 $file ${OUTPREFIX}5_cluster_aligned/$bname.sam $PROTUMIREF > ${OUTPREFIX}5_racon/$bname.fasta
      # polish with medaka
      ${medaka_consensusLOC} -i $file -d ${OUTPREFIX}5_racon/$bname.fasta -o ${OUTPREFIX}5_medaka/$bname -t 1 -m r941_min_hac_g507 
      # align clusters to polished references
      ${minimap2LOC} -ax map-ont -t 1 --secondary=no ${OUTPREFIX}5_medaka/$bname/consensus.fasta $file | ${samtoolsLOC} view -hF 2048 > ${OUTPREFIX}5_cluster_polished_aligned/$bname.sam
}
export -f dowork

# get cluster files
files=$(find ${OUTPREFIX}4_cluster_fastq/ -name "*.fastq")

# parallel processing of function
${parallelLOC} -k -j $THREADS dowork {} ${PROTUMIREF} ${OUTPREFIX} ${MODEL} ${minimap2LOC} ${samtoolsLOC} ${medaka_consensusLOC} ${raconLOC} ::: $files


# write table header
echo -e "Cluster\tID\tProtein_score" > ${OUTPREFIX}5_cluster_polished_aligned/cluster_eval.tsv

# loop through all clusters and collect sequences and alignment scores
for file in $files
do
      bname=$(basename -s .fastq $file)
      # polished sequences
      echo ">$bname" >> ${OUTPREFIX}polished_sequences.fasta
      sed -n 2p ${OUTPREFIX}5_medaka/$bname/consensus.fasta >> ${OUTPREFIX}polished_sequences.fasta
      # alignment scores
      ${samtoolsLOC} view ${OUTPREFIX}5_cluster_polished_aligned/$bname.sam | awk -v bname=$bname -F '\t' '{ OFS="\t" } {print bname,$1,$14}' >> ${OUTPREFIX}5_cluster_polished_aligned/cluster_eval.tsv
done

# align all polished sequences to the references
${minimap2LOC} -a --secondary=no -O 8,28 $PROTUMIREF ${OUTPREFIX}polished_sequences.fasta | ${samtoolsLOC} view -hbF 2048 | ${samtoolsLOC} sort > ${OUTPREFIX}polished_sequences.bam
${samtoolsLOC} index ${OUTPREFIX}polished_sequences.bam
