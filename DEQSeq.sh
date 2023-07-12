#!/bin/bash
# DEQSeq pipeline; cluster and polish nanopore reads with UMI while counting DNA editing on the target site of the same read
# Lukas T. Schmitt, Buchholz Lab, Medical Faculty, TU Dresden

######## argument passing ########

function usage () {
    cat >&2 <<EOF
  USAGE:  [options]
	-h  Print the usage info.
	-y  <YAML config file>	: Path to the YAML config file. Required.
	-s  <start from nr>	: Step to start from, integer. Default: 0
	-d  <DEQSeq directory>	: Directory containing DEQSeq scripts. Default: path to this script.
EOF
}

# Get default script directory based on link location
SCRIPTDIR=$(dirname $(readlink -f $0))

# start from step 0 if not told otherwise
STARTFROM=0

# Get input variables #
while getopts "hy:d:s:" options; do
      case ${options} in
      h ) usage
	  exit 1;;
      y ) yaml=${OPTARG};;
      s ) STARTFROM=${OPTARG};;
      d ) SCRIPTDIR=${OPTARG};;
      esac
done

# Print usage if no variables are entered
if [[ ${OPTIND} -eq 1 ]] ; then
    usage
    exit 1
fi

######## get application locations ########

#function read_check_loc () {
#      $1 = $yaml_name
#      $2 = $variable
#      if [ -z $variable ]; then
#            echo "Argument ${yaml_name} in yaml not set!" > stderr
#	    exit 1
#      fi
#      # check if file/program exists
#      # output value for assignment
#}

if grep -q 'samtools_location:' ${yaml} ; then
    samtoolsLOC=$(grep 'samtools_location:' ${yaml} | awk '{print $2}')
else
    samtoolsLOC=samtools
fi

if grep -q 'Rscript_location:' ${yaml} ; then
    RLOC=$(grep 'Rscript_location:' ${yaml} | awk '{print $2}')
else
    RLOC=Rscript
fi

if grep -q 'minimap2_location:' ${yaml} ; then
    minimap2LOC=$(grep 'minimap2_location:' ${yaml} | awk '{print $2}')
else
    minimap2LOC=minimap2
fi

if grep -q 'parallel_location:' ${yaml} ; then
    parallelLOC=$(grep 'parallel_location:' ${yaml} | awk '{print $2}')
else
    parallelLOC=parallel
fi

if grep -q 'racon_location:' ${yaml} ; then
    raconLOC=$(grep 'racon_location:' ${yaml} | awk '{print $2}')
else
    raconLOC=racon
fi

if grep -q 'medaka_consensus_location:' ${yaml} ; then
    medaka_consensusLOC=$(grep 'medaka_consensus_location:' ${yaml} | awk '{print $2}')
else
    medaka_consensusLOC=medaka_consensus
fi

if grep -q 'vsearch_location:' ${yaml} ; then
    vsearchLOC=$(grep 'vsearch_location:' ${yaml} | awk '{print $2}')
else
    vsearchLOC=vsearch
fi

######## get script variables ########

OUTPREFIX=$(grep 'output_prefix:' ${yaml} | awk '{print $2}')
READS=$(grep 'fastq_reads:' ${yaml} | awk '{print $2}')
THREADS=$(grep 'num_threads:' ${yaml} | awk '{print $2}')

PROTUMIREF=$(grep 'protein_umi_fasta:' ${yaml} | awk '{print $2}')
UMILOC=$(grep 'umi_locations:' ${yaml} | awk '{print $2}')
PROTLOC=$(grep 'protein_locations:' ${yaml} | awk '{print $2}')
TSREF=$(grep 'target_sites_fasta:' ${yaml} | awk '{print $2}')
TSLOC=$(grep 'target_site_locations:' ${yaml} | awk '{print $2}')

MINCLUSTSIZE=$(grep 'minimum_cluster_size:' ${yaml} | awk '{print $2}')

UMIIDENT1=$(grep 'umi_cluster_identity1:' ${yaml} | awk '{print $2}')
UMIIDENT2=$(grep 'umi_cluster_identity2:' ${yaml} | awk '{print $2}')

MODEL=$(grep 'medaka_model:' ${yaml} | awk '{print $2}')

CLEANUP=$(grep 'cleanup_lvl:' ${yaml} | awk '{print $2}')

######## run scripts ########

# make reference bed files for alignment filtering
if [[ ${STARTFROM} -lt 1 ]]; then
      ${RLOC} ${SCRIPTDIR}/0_make-references.R ${OUTPREFIX} ${PROTUMIREF} ${TSREF} || exit 1
fi

# filter reads and align to Protein UMI reference
if [[ ${STARTFROM} -lt 2 ]]; then
      echo "####### 1. Preprocessing #######"
      bash ${SCRIPTDIR}/1_preprocessing.sh ${OUTPREFIX} ${PROTUMIREF} ${READS} ${THREADS} ${minimap2LOC} ${samtoolsLOC} || exit 1
fi

# extract UMIs
if [[ ${STARTFROM} -lt 3 ]]; then
      echo "####### 2. UMI extraction #######"
      ${RLOC} ${SCRIPTDIR}/2_UMI-extraction.R ${OUTPREFIX} ${UMILOC} || exit 1
fi

# cluster UMIs
if [[ ${STARTFROM} -lt 4 ]]; then
      echo "####### 3. Clustering #######"
      bash ${SCRIPTDIR}/3_cluster.sh ${OUTPREFIX} ${UMIIDENT1} ${UMIIDENT2} ${THREADS} ${vsearchLOC} || exit 1
fi

# split reads according to clusters
if [[ ${STARTFROM} -lt 5 ]]; then
      echo "####### 4. Split clusters #######"
      ${RLOC} ${SCRIPTDIR}/4_split-clusters.R ${OUTPREFIX} ${MINCLUSTSIZE} ${UMIIDENT2} || exit 1
fi

# make accurate dna sequences
if [[ ${STARTFROM} -lt 6 ]]; then
      echo "####### 5. Polish sequences #######"
      bash ${SCRIPTDIR}/5_polish.sh ${OUTPREFIX} ${PROTUMIREF} ${THREADS} ${MODEL} ${minimap2LOC} ${samtoolsLOC} ${medaka_consensusLOC} ${raconLOC} ${parallelLOC} || exit 1
fi

# align to possible editing outcomes
if [[ ${STARTFROM} -lt 7 ]]; then
      echo "####### 6. Get edits #######"
      bash ${SCRIPTDIR}/6a_edit_match.sh ${OUTPREFIX} ${TSREF} ${THREADS} ${minimap2LOC} ${samtoolsLOC} ${parallelLOC} || exit 1
      if [[ -n ${TSLOC} ]]; then
	    echo "Extracting edit region..."
	    ${RLOC} ${SCRIPTDIR}/6b_edit_extract.R ${OUTPREFIX} ${TSREF} ${TSLOC} || exit 1
      fi
fi

# evaluate, extract protein, count edits and combine data
if [[ ${STARTFROM} -lt 8 ]]; then
      echo "####### 7. Summarise #######"
      ${RLOC} ${SCRIPTDIR}/7_extract_count.R ${OUTPREFIX} ${PROTUMIREF} ${UMILOC} ${PROTLOC} ${MINCLUSTSIZE} || exit 1
fi

if [[ ${STARTFROM} -lt 9 ]]; then
      echo "####### 8. Cleanup #######"
      # level 1 and level 2 delete, 
      # 1: all folders except 5_cluser_aligned and 6_ts_aligned
      # 2: also delete these two folders as well as the polished_sequences.bam + index
      if [[ ${CLEANUP} -gt 0 ]]; then
	    rm -rv ${OUTPREFIX}{0,1,2,4}_* ${OUTPREFIX}5_cluster_polished_aligned ${OUTPREFIX}5_medaka ${OUTPREFIX}5_racon
      fi
      if [[ ${CLEANUP} -gt 1 ]]; then
	    rm -rv ${OUTPREFIX}polished_sequences.bam* ${OUTPREFIX}5_cluster_aligned ${OUTPREFIX}6_ts_aligned
      fi
fi

echo "####### Done! #######"
