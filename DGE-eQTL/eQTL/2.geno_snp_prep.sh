#!/bin/bash
set -ux

# Aims:
# Gather snp and genotyping data and make meQTL files
# run script as source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh

# VAST vars
plink_data=vast
dir=~/GWAS_22/gwas_final/merge/typhoid/vast
plink_dir=${dir}
outdir=~/GWAS_22/gwas_final/eQTL/vast/nov8
study=vast_
data=V1WARS_iga

# plink_data=typhoid2.IBD
# plink_dir=~/GWAS_22/gwas_final/merge/typhoid/QC
# outdir=~/GWAS_22/gwas_final/eQTL/vast_tyg/Nov/combat/redo
# data=TD_just_T_snpG_
# study=vast_tyg

# same as outdir in Rscript :: EQTL = ~/GWAS_22/gwas_final/eQTL # Or adnob/antibody measure
cd ${outdir}
log=${outdir}/log

# Make snp keep file from gwas data (not needed for eqtl)
# cut -f3 tophits_${data}.txt > ${study}${data}_snp_keep.txt

# Filter plink data for study individuals, and snps for eqtl analysis
plink2 --bfile ${plink_dir}/${plink_data} \
--keep ${data}_keep.txt \
--extract ${data}_snp_keep.txt \
--make-bed \
--out ${data}

# # Filter for high LD SNPs
# 
# plink2 --bfile ${data} \
# --indep-pairwise 50 5 0.95 \
# --out ${data}.LD \
# --bad-ld
# 
# # Prune variants
# plink2 \
# --bfile ${data} \
# --allow-extra-chr \
# --extract ${data}.LD.prune.in \
# --make-bed \
# --out ${data}.LD

#Recode SNPs as 1's and 0's
# Transpose table so SNPs are rows, participants are columns
plink --bfile ${data} \
--allow-no-sex \
--recode oxford \
--out ${data}

plink \
--data ${data} \
--allow-no-sex \
--recode A-transpose \
--out ${data}

#  CHR	SNP	(C)M	POS	COUNTED	ALT
# Cut relevant data from new file
cat ${data}.traw | cut -f2,7- > ${data}_geno.txt
awk '{print $2, $1, $4}' ${data}.traw > ${data}_snp_loc.txt

rm *.nosex
mv *.log ${log}