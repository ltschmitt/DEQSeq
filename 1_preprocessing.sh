#!/bin/bash

OUTPREFIX=$1
PROTUMIREF=$2 
READS=$3 
MINLEN=$4 
MINQUAL=$5 
THREADS=$6
filtlongLOC=$7
minimap2LOC=$8
samtoolsLOC=$9

mkdir -p ${OUTPREFIX}1_preprocessing

# filter reads by size
${filtlongLOC} --min_length $MINLEN --min_mean_q $MINQUAL $READS > ${OUTPREFIX}1_preprocessing/filtered.fastq
# align reads to references and filter for reads spanning most of the reference
${minimap2LOC} -t $THREADS -ax map-ont --secondary no $PROTUMIREF ${OUTPREFIX}1_preprocessing/filtered.fastq | ${samtoolsLOC} sort | ${samtoolsLOC} view -h -F 2048 -L ${OUTPREFIX}0_references/protumi1.bed | ${samtoolsLOC} view -hb -L ${OUTPREFIX}0_references/protumi2.bed > ${OUTPREFIX}1_preprocessing/regionfiltered.bam
${samtoolsLOC} index ${OUTPREFIX}1_preprocessing/regionfiltered.bam







