#!/bin/bash
set -ux

# Filter topmed imputation data, index, concatanate and convert to plink

# Steps
# 1) Download data (unzip topmed file with password)
# 2) Filter for only high qual SNPs R2 > 0.8, filter for multiallelic snps
# 3) Convert to bcf and index
# 4) Check file against index
# 5) annotate CHR POS REF
# 6) merge

# Command info at https://samtools.github.io/bcftools/bcftools.html
# Run as source ~/GWAS_22/Enteric_GWAS/Post-Impute/QC/impute_filter.sh

# Directory and Variables -----------------------------------------------
DIR=~/GWAS_22/new_gwas/Post-imp/VAST
OUT_DIR=~/GWAS_22/gwas_final/VAST
REF_DIR=GWAS/TOPMED_impute
FILENAME=VAST

# Post-impute filtering of vcf.gz files --------------------------------
for chr_n in {1..22}; do 
bcftools view -Ou -i 'R2>0.8' ${DIR}/chr${chr_n}.dose.vcf.gz |  #R2 filter 
bcftools view -Ou -m2 -M2 -v snps |  # Multi-allelic filter
bcftools norm -Ou -m -any |  #Index
bcftools norm -Ou -f ${REF_DIR}/GRCh38.fa |
bcftools annotate -Ou -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ou |  #New Chr name
bcftools norm --rm-dup all -Ob ${OUT_DIR}/chr${chr_n}.index.R2_08.bcf
echo File indexing and filtering DONE

bcftools concat --naieve ${OUT_DIR}/chr{1..22}.index.R2_08.bcf -Ob -o ${FILENAME}_concatR2.bcf
echo bcf Chr concat DONE

plink2 --bcf ${FILENAME}_concat.bcf \
--keep-allele-order \
--double-id \
--biallelic-only \
--make-bed \
--allow-extra-chr 0 \

echo Plink conversion DONE