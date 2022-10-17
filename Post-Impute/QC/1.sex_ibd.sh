#!/bin/bash

# IBD script checking
# source ~/GWAS_22/Enteric_GWAS/Post-Impute/QC/1.sex_ibd.sh
set -ux

# Variables and Directories ---------------------------------------------------------

dir=~/GWAS_22/gwas_final/merge
data=typhoid

cd ${dir}
mkdir QC
qcdir=${dir}/QC
mkdir log
log=${dir}/log

cd ${qcdir}

# Filter failed sex check samples ---------------------------------------------------

plink2 \
--bfile ${dir}/${data} \
--remove ${dir}/fail_sex.txt \
--make-bed \
--out ${dir}/${data}

# Calculate IBD in Plink ------------------------------------------------------------

plink \
--bfile ${dir}/${data} \
--genome \
--out ${qcdir}/${data}

echo "IBD genome - Done"

# Remove one sample from each pair with pi-hat (% IBD) above threshold (0.1875):
awk '$10 >= 0.1875 {print $2, $4, $10}' ${data}.genome > pair_list.txt
awk '$10 >= 0.1875 {print $1, $2}' ${data}.genome | uniq > ${data}.outliers.txt 
wc -l ${data}.outliers.txt

echo "Outlier list - Done"

# Filter IBD outliers in PLINK -------------------------------------------------------
plink2 \
--bfile ${dir}/${data} \
--remove ${data}.outliers.txt \
--make-bed \
--out ${qcdir}/${data}.IBD

echo "Related Samples Removed"