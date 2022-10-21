#!/bin/bash

set -ux

# Aims:
# Gather snp and genotyping data and make meQTL files
# run script as source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh
plink_data=typhoid2.IBD
plink_dir=~/GWAS_22/gwas_final/merge/typhoid/QC

outdir=~/GWAS_22/gwas_final/eQTL # same as outdir in Rscript
time_point=12h
study=T1T2_

cd ${outdir}

# Filter plink data for study individuals, and snps for eqtl analysis
plink2 --bfile ${plink_dir}/${plink_data} \
--keep ${study}${time_point}_keep.txt \
--extract ${study}${time_point}_snp_keep.txt \
--make-bed \
--out ${study}${time_point}

#Recode SNPs as 1's and 0's
# Transpose table so SNPs are rows, participants are columns
plink --bfile ${study}${time_point} \
--allow-no-sex \
--recode oxford \
--out ${study}${time_point}

plink --data ${study}${time_point} \
--allow-no-sex \
--recode A-transpose \
--out ${study}${time_point}

#  CHR	SNP	(C)M	POS	COUNTED	ALT
# Cut relevant data from new file
cat ${study}${time_point}.traw | cut -f2,7- > ${study}${time_point}_geno.txt
awk '{print $2, $1, $4}' ${study}${time_point}.traw > ${study}${time_point}_snp_loc.txt

rm *.nosex