#!/bin/bash
#
# Generate knockoff negative-controls
#
# Authors: Matteo Sesia
# Date:    07/21/2020

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR
mkdir -p $TMP_DIR"/knockoffs"

# Storage of output files
OUT_DIR="../results"
mkdir -p $OUT_DIR

# List of chromosomes
# SA hard coding chr22 and commenting out examples 
CHR=22 # im not using this anywhere but down the road would put ${CHR} in spot where i have chr22
# CHR_MIN=21
# CHR_MAX=22

# Path to snpknock2 executable built as described above
SNPKNOCK2="../snpknock2/bin/snpknock2"

# Which operations should we perform?
FLAG_GENERATE_KNOCKOFFS=1

######################
# Generate knockoffs #
######################

if [[ $FLAG_GENERATE_KNOCKOFFS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Generating knockoffs"
  echo "----------------------------------------------------------------------------------------------------"

  # Run snpknock2
# SA i removed the ibd line - it was under the --part
# --ibd "../data/ibd/example_chr{"$CHR_MIN":"$CHR_MAX"}.txt" \


  $SNPKNOCK2 \
    --haps ../data/haplotypes/unicorn/chr22_Ana_haps \
    --keep "../data/qc/chr22_Ana_samples_qc.txt" \
    --extract "../data/qc/qc_chr22_Ana_only_colons.txt" \
    --map "../data/maps/chr22_Ana_with_rate_fixed.b37.gmap" \
    --part "../tmp/partitions/chr22_Ana_only_colons.txt" \
    --K 10 \
    --cluster_size_min 1000 \
    --cluster_size_max 10000 \
    --hmm-rho 1 \
    --hmm-lambda 1e-3 \
    --windows 0 \
    --n_threads 4 \
    --seed 2020 \
    --compute-references \
    --generate-knockoffs \
    --out $TMP_DIR"/knockoffs/chr22_Ana/knockoffs_chr22_Ana"

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping generation of knockoffs"
  echo "----------------------------------------------------------------------------------------------------"
fi