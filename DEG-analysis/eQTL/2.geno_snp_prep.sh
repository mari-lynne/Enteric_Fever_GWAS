#!/bin/bash

set -ux

# Aims:
# Gather snp and genotyping data and make meQTL files
# run script as source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh
name=T1T2
dir=~/GWAS_22/gwas_final/eQTL

cd ${dir}
#Recode SNPs as 1's and 0's
# Transpose table so SNPs are rows, participants are columns
plink --bfile ${name} \
--allow-no-sex \
--recode oxford \
--out ${name}

plink --data ${name} \
--allow-no-sex \
--recode A-transpose \
--out ${name}

#  CHR	SNP	(C)M	POS	COUNTED	ALT
# Cut relevant data from new file
cat ${name}.traw | cut -f2,7- > ${name}_geno_file.txt
awk '{print $2, $1, $4}' ${name}.traw > ${name}_snp_loc.txt

rm *.nosex