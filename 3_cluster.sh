#!/bin/bash
OUTPREFIX=$1
UMIIDENT1=$2
UMIIDENT2=$3
THREADS=$4
vsearchLOC=$5

${vsearchLOC} --threads $THREADS --id $UMIIDENT1 --cluster_size ${OUTPREFIX}2_3_vsearch/umis.fasta --consout ${OUTPREFIX}2_3_vsearch/consout.fa --clusterout_id --clusterout_sort --uc ${OUTPREFIX}2_3_vsearch/uc.tab

# only run second clustering if UMIIDENT2 is a value between 0-1
if echo $UMIIDENT2 | awk '($1 > 1 || $1 < 0) { exit 1 }'; then
      ${vsearchLOC} --threads $THREADS --id $UMIIDENT2 --cluster_size ${OUTPREFIX}2_3_vsearch/consout.fa --clusterout_id --clusterout_sort --uc ${OUTPREFIX}2_3_vsearch/uc2.tab
fi


