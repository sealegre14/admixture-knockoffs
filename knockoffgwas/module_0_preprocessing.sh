erm why did the 0 id_id print to my terminal when i ran this 

#!/bin/bash
#
# Preprocess VCF data --> take in VCF --> assign unique IDs 
# --> PLINk for BED files --> PLINK for HAPS files --> legened and qc file creation
# 
# Author: Sarah Alegre
# Date:    02/18/2026


set -e  # stop if anything breaks

# instructions for expected input - expects atleast 1 VCF input 
# ./module_0_preprocessing path/to/VCF/chr22_Ana.vcf
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 file1.vcf [file2.vcf ...]"
  exit 1
fi

# allow the user to name current batch of VCFs (first argument)
BATCH_NAME="$1"
shift  # remove batch name from arguments, now "$@" contains only VCFs

OUT_DIR="../data/clean_vcf"  # output directory
mkdir -p $OUT_DIR

QC_DIR="../data/qc/$BATCH_NAME"  # QC directory for this batch
mkdir -p "$QC_DIR"
echo "QC files will be written to $QC_DIR"

# clean up old 
> "$QC_DIR/variants_qc.txt"   # truncate or create empty
> "$QC_DIR/samples_qc.txt"    # truncate or create empty


MAP_DIR="../data/maps"       
UPDATED_MAP_DIR="$MAP_DIR/$BATCH_NAME"  # map directory for this batch
mkdir -p "$UPDATED_MAP_DIR"

# loop for cleaning VCF and adding unique IDs
for VCF in "$@"; do

  echo "Processing $VCF..."

  BASENAME=$(basename "$VCF" .vcf)
  OUT_FILE="$OUT_DIR/${BASENAME}_clean.vcf"

  # Step 1: Clean VCF 
  # rewrite ID name to be CHROM:POS:REF:ALT
  bcftools annotate \
    --set-id '%CHROM:%POS:%REF:%ALT' \
    "$VCF" \
    -Ov -o "$OUT_FILE"

  echo "Created $OUT_FILE"

  # Step 2: Make output directories for VCF
  GENO_DIR="../data/genotypes/${BASENAME}"  # genotypes directory and subdirectory
  mkdir -p "$GENO_DIR"
  HAP_DIR="../data/haplotypes/${BASENAME}"  # haplotypes directory and subdirectory 
  mkdir -p "$HAP_DIR"

  # Step 3: PLINK → BED
  plink \
    --vcf "$OUT_FILE" \
    --make-bed \
    --out "$GENO_DIR/${BASENAME}" \
    --allow-extra-chr

  # Step 4: PLINK2 → HAPS
  plink2 \
    --vcf "$OUT_FILE" \
    --export haps \
    --out "$HAP_DIR/${BASENAME}"

  # Step 5: generate qc files - variants, samples, and qc per vcf 
  GENO_BIM="$GENO_DIR/${BASENAME}.bim"  # where we just generated the BIM files

  # Per-VCF variant QC
  awk '{print $2}' "$GENO_BIM" > "$QC_DIR/qc_${BASENAME}.txt"

  # Append batch-wide variant IDs
  awk '{print $2}' "$GENO_BIM" >> "$QC_DIR/variants_qc.txt"

  # Append batch-wide sample IDs (matching .sample)
  awk 'NR>2 {print $1, $2}' "$HAP_DIR/${BASENAME}.sample" >> "$QC_DIR/samples_qc.txt"

  echo "Finished $BASENAME"

done

# Step 6: calculate recombination rate in gmap files, and update format 
# Loop over all chromosomes for which VCFs were processed
for VCF in "$@"; do
    BASENAME=$(basename "$VCF" .vcf)
    
    # Extract chromosome name from VCF basename (assumes format chrXX_*)
    CHR=$(echo "$BASENAME" | cut -d'_' -f1)
    
    ORIG_MAP="$MAP_DIR/${CHR}.b37.gmap"
    UPDATED_MAP="$UPDATED_MAP_DIR/${CHR}.b37.gmap.clean.txt"
    
    if [ ! -f "$ORIG_MAP" ]; then
        echo "WARNING: Map file not found: $ORIG_MAP"
        continue
    fi
    
    echo "Updating map for $CHR..."
    
    awk 'NR==1 {prev_cM=$3; prev_BP=$1; rate=0; print $2,$1,rate,$3; next}
    {
      delta_cM = $3 - prev_cM
      delta_BP = $1 - prev_BP
      if(delta_BP==0){rate=0}else{rate = delta_cM/(delta_BP/1e6)}
      print $2,$1,rate,$3
      prev_cM=$3; prev_BP=$1
    }' "$ORIG_MAP" > "$UPDATED_MAP"

    
    echo "Updated map written to $UPDATED_MAP"
done



# Function to check duplicates in a file
check_duplicates() {
    local FILE="$1"
    local DESC="$2"
    local DUPES
    DUPES=$(sort "$FILE" | uniq -d)
    if [ -n "$DUPES" ]; then
        echo "WARNING: duplicate $DESC found in $FILE!"
        echo "$DUPES"
    fi
}

# Check duplicates for batch-wide QC files
check_duplicates "$QC_DIR/variants_qc.txt" "variant IDs"
check_duplicates "$QC_DIR/samples_qc.txt" "sample IDs"


echo "Step 1 complete."