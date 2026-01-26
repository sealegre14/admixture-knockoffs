#!/usr/bin/env Rscript

# ===============================
# partition_binned.R
# ===============================
# Partition SNPs into fixed genetic distance bins (5 cM and 10 cM)
# Replacement for the original partition.R
# Compatible with Module 1 shell script and downstream modules
# ===============================

suppressMessages(library(tidyverse))

# -------------------------------
# Read input arguments
# -------------------------------
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4){
  stop("Usage: Rscript partition_binned.R <map.file> <bim.file> <qc.file> <out.file>")
}

map.file  <- args[1]
bim.file  <- args[2]
qc.file   <- args[3]
out.file  <- args[4]

cat(sprintf("Map file: %s\nBIM file: %s\nQC file: %s\nOutput file: %s\n\n",
            map.file, bim.file, qc.file, out.file))

# -------------------------------
# Load BIM file (SNP positions)
# -------------------------------
Variants <- read_delim(
  bim.file,
  delim = "\t",
  col_names = c("CHR", "SNP", "X", "BP", "A1", "A2"),
  col_types = cols()
) %>%
  select(CHR, SNP, BP)

# -------------------------------
# Load QC variants
# -------------------------------
variants.qc <- read_delim(
  qc.file,
  delim = " ",
  col_names = "SNP",
  col_types = cols()
)

# Keep only QC-passed SNPs
Variants <- Variants %>% inner_join(variants.qc, by = "SNP")

# -------------------------------
# Load genetic map
# -------------------------------
Map <- read_delim(
  map.file,
  delim = "\t",
  col_types = cols()
) %>%
  transmute(
    CHR = as.integer(chr),
    BP  = as.integer(pos),
    cM  = as.numeric(cM)
  )

# -------------------------------
# Cross-reference SNPs with map
# -------------------------------
Map <- Map %>% inner_join(Variants, by = c("CHR", "BP"))

# -------------------------------
# Assign bins for multiple resolutions
# -------------------------------
# 2 resolutions: 5 cM and 10 cM
resolutions <- c(5, 10)

# Compute bin IDs for each SNP
Partitions <- Map %>% select(SNP, cM)
for(i in seq_along(resolutions)){
  bin_size <- resolutions[i]
  col_name <- paste0("res_", i)
  # floor division, plus 1 so first bin = 1
  Partitions[[col_name]] <- floor(Partitions$cM / bin_size) + 1
}

# -------------------------------
# Save output
# -------------------------------
Partitions %>% select(-cM) %>%
  write_delim(out.file, delim = " ")

cat(sprintf("Partitioning complete. Output written to: %s\n", out.file))
