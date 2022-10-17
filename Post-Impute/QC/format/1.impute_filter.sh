#!/bin/bash
set -ux

# Filter topmed imputation data, index, concatanate and convert to plink

# Steps
# 1) Download data (unzip topmed file with password)
# 2) Filter for only high qual SNPs R2 > 0.8, filter for multiallelic snps
# 3) Convert to bcf and index
# 4) Check file against index
# 5) annotate CHR POS REF

# Command info at https://samtools.github.io/bcftools/bcftools.html
# Run as source ~/GWAS_22/Enteric_GWAS/Post-Impute/QC/format/2.impute_filter.sh

# Directory and Variables -----------------------------------------------
DIR=~/GWAS_22/new_gwas/Post-imp/P1T1
OUT_DIR=~/GWAS_22/gwas_final/P1T1
REF_DIR=~/GWAS/TOPMED_impute
FILENAME=P1T1

# Post-impute filtering of vcf.gz files --------------------------------
for chr_n in {1..22}; do 

 echo "Processing Chromosome ${chr_n}"
 echo "======================================"

bcftools view -Ou -i 'R2>0.8' ${DIR}/chr${chr_n}.dose.vcf.gz \
| bcftools view -Ou -m2 -M2 -v snps \
| bcftools norm -Ou -m -any \
| bcftools norm -Ou -f ${REF_DIR}/GRCh38.fa \
| bcftools annotate -Ou -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ou \
| bcftools norm --rm-dup all -Ob > ${OUT_DIR}/chr${chr_n}.index.R2_08.bcf; done
echo "Chromosome${chr_n} index and filtering DONE"

# Update headers
for chr_n in {1..22}; do
sed -i "s/contig=<ID=${chr_n}/contig=<ID=chr${chr_n}/g" ${OUT_DIR}/chr${chr_n}.index.R2_08.bcf; done
echo "Chr${chr_n} header updating DONE"

# Concat
bcftools concat ${OUT_DIR}/chr{1..22}.index.R2_08.bcf -Ob -o ${OUT_DIR}/${FILENAME}_concatR2.bcf
echo "all chromosomes bcf-concat DONE"