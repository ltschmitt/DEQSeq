#!/bin/bash
OUTPREFIX=$1
TSREF=$2
THREADS=$3
minimap2LOC=$4
samtoolsLOC=$5
parallelLOC=$6

mkdir ${OUTPREFIX}6_ts_aligned

# function for parallel mapping
dowork() {
      OUTPREFIX=$1
      TSREF=$2
      file=$3
      minimap2LOC=$4
      samtoolsLOC=$5
      bname=$(basename -s .fastq $file)
      
      ${minimap2LOC} -t 1 --secondary=no -ax map-ont $TSREF $file | ${samtoolsLOC} view -hF 2048 -L ${OUTPREFIX}0_references/ts1.bed | ${samtoolsLOC} view -hbL ${OUTPREFIX}0_references/ts2.bed | ${samtoolsLOC} sort > ${OUTPREFIX}6_ts_aligned/$bname.bam
      ${samtoolsLOC} index ${OUTPREFIX}6_ts_aligned/$bname.bam
}
export -f dowork

# get cluster files
files=$(find ${OUTPREFIX}4_cluster_fastq/ -name "*.fastq")

# parallel processing of function
${parallelLOC} -k -j $THREADS dowork $1 $2 {} $4 $5 ::: $files


# make outputfile
echo -e "Cluster\tID\tTS\tTS_score" > ${OUTPREFIX}6_ts_aligned/ts_editing.tsv

# gather alignment information
for file in `find ${OUTPREFIX}4_cluster_fastq -name "*.fastq"`
do
      bname=$(basename -s .fastq $file) 
      ${samtoolsLOC} view ${OUTPREFIX}6_ts_aligned/$bname.bam | awk -v bname=$bname -F '\t' '{ OFS="\t"} {print bname,$1,$3,$14}'
done >> ${OUTPREFIX}6_ts_aligned/ts_editing.tsv
