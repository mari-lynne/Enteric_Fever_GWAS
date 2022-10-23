#!/bin/bash

# Remove missex/missing data samples
# Calcualte IBD
# Remove cryptic related samples

# source ~/GWAS_22/Enteric_GWAS/Post-Impute/QC/1.sex_ibd.sh
set -ux

# Variables and Directories ---------------------------------------------------------

dir=~/GWAS_22/gwas_final/merge
plink_data=merge_rn
data=enteric

cd ${dir}
#mkdir QC
qcdir=${dir}/enteric/QC
mkdir log
log=${qcdir}/log

cd ${qcdir}

# Filter samples ---------------------------------------------------

# Miss-sex and missing data samples filtered in pheno_covar.R #fail_sex.txt

plink2 \
--bfile ${dir}/${plink_data} \
--keep ${dir}/keep_ent.txt \
--make-bed \
--out ${qcdir}/${data}

# Calculate IBD in Plink ------------------------------------------------------------

plink \
--bfile ${data} \
--genome \
--out ${data}

echo "IBD genome - Done"

# Remove one sample from each pair with pi-hat (% IBD) above threshold (0.1875):
awk '$10 >= 0.1875 {print $2, $4, $10}' ${data}.genome > pair_list.txt
awk '$10 >= 0.1875 {print $1, $2}' ${data}.genome | uniq > ${data}.outliers.txt 
wc -l ${data}.outliers.txt

echo "Outlier list - Done"

# Filter IBD outliers in PLINK -------------------------------------------------------
plink2 \
--bfile ${data} \
--remove ${data}.outliers.txt \
--make-bed \
--out ${data}.IBD

echo "Related Samples Removed"

mv *.log ${log}
rm -f *.nosex