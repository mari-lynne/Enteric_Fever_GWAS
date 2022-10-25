#!/bin/bash
set -ux

# Aims:
# Gather snp and genotyping data and make meQTL files
# run script as source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh
plink_data=typhoid2.IBD
dir=~/GWAS_22/gwas_final/merge/typhoid
plink_dir=${dir}/QC

#~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus/adnob
# data renamed to data
data=0.5
study=vast_tyg_
outdir=~/GWAS_22/gwas_final/eQTL
# same as outdir in Rscript :: EQTL = ~/GWAS_22/gwas_final/eQTL
 # Or adnob/antibody measure

cd ${outdir}
# Make snp keep file from gwas data (not needed for eqtl)
# cut -f3 tophits_${data}.txt > ${study}${data}_snp_keep.txt

# Filter plink data for study individuals, and snps for eqtl analysis
plink2 --bfile ${plink_dir}/${plink_data} \
--keep ${study}${data}_keep.txt \
--extract ${study}${data}_snp_keep.txt \
--make-bed \
--out ${study}${data}

#Recode SNPs as 1's and 0's
# Transpose table so SNPs are rows, participants are columns
plink --bfile ${study}${data} \
--allow-no-sex \
--recode oxford \
--out ${study}${data}

plink --data ${study}${data} \
--allow-no-sex \
--recode A-transpose \
--out ${study}${data}

#  CHR	SNP	(C)M	POS	COUNTED	ALT
# Cut relevant data from new file
cat ${study}${data}.traw | cut -f2,7- > ${study}${data}_geno.txt
awk '{print $2, $1, $4}' ${study}${data}.traw > ${study}${data}_snp_loc.txt

rm *.nosex
mv *log log