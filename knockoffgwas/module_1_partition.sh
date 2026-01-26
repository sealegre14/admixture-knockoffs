#!/bin/bash
#
# Partition the variants through adjacent clustering
#
# Authors: Matteo Sesia
# Date:    07/21/2020

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR
mkdir -p $TMP_DIR"/partitions"

# List of chromosomes
# CHR_LIST=$(seq 21 22)
CHR_LIST="22"

# List of resolutions in cM (note: this is also fixed inside the R script)
RESOLUTION_LIST=("0" "0.01" "0.05" "0.1" "0.2" "0.5" "1")


# Utility scripts
#SA - changing to alternate partition_binned.R script for binning instead of clustering 
PARTITION_VARIANTS="Rscript --vanilla utils/partition.R"
# PARTITION_VARIANTS="Rscript --vanilla utils/partition_binned.R"
# SA now trying python script instead 
# PARTITION_VARIANTS="utils/partition_binned_numpy.py"

# Which operations should we perform?
FLAG_PARTITION=1

################
# Partitioning #
################

if [[ $FLAG_PARTITION == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Partitioning variants"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do

    echo ""
    echo "Processing chromosome "$CHR" ..."
    echo ""

    # Input genotype files (PLINK format)
    # SA - hard coding to chr22_Ana
    # GENO_BIM="../data/genotypes/example_chr"$CHR".bim"
    GENO_BIM="../data/genotypes/chr22_Ana.bim"

    # List of variants that passed QC
    # SA giving fake qc file for now 
    # QC_VARIANTS="../data/qc/variants_qc.txt"
    QC_VARIANTS="../data/qc/chr22_bad_id_qc_mod1.txt"

    # SA - adding a block to generate qc variant file if it doesnt exist - one long column with matching IDs
    if [ ! -s "$QC_VARIANTS" ]; then
      echo "QC file missing or empty — generating from BIM..."
      mkdir -p ../data/qc
      awk '{print $2}' "$GENO_BIM" > "$QC_VARIANTS"
    fi

    # Genetic map file
    # SA changing this to map to the .gmap version which i downloaded from github cole gave me 
    #GEN_MAP="../data/maps/genetic_map_chr"$CHR".txt"
    GEN_MAP="../data/maps/chr22.b37.gmap"

    # Basename for output dendrogram file
    # SA - hard coding to chr22_Ana
    # OUT_FILE=$TMP_DIR"/partitions/example_chr"$CHR".txt"
    OUT_FILE=$TMP_DIR"/partitions/chr22_Ana_PLEASE.txt"

    # Partition the variants at different resolutions
    $PARTITION_VARIANTS $GEN_MAP $GENO_BIM $QC_VARIANTS $OUT_FILE
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping variant partitioning"
  echo "----------------------------------------------------------------------------------------------------"
fi
