#!/bin/bash
set -ux
# Subset data and rerun LD/PCA and make covar file

# source ~/GWAS_22/Enteric_GWAS/Post-Impute/QC/study_pca.sh

# Variables and Directories -----------------------------------------------------------
#plink_data=typhoid2
data=para

dir=~/GWAS_22/gwas_final/merge
qcdir=${dir}/${data}/QC

cd ${qcdir}

mkdir pca
pcadir=${qcdir}/pca

# Filter data for study samples
# plink2 \
# --bfile ${dir}/${plink_data} \
# --keep ${qcdir}/${data}_keep.txt \
# --make-bed \
# --out ${qcdir}/${data}

# Calculate LD

plink2 \
--bfile ${qcdir}/${data}.IBD \
--indep-pairwise 50 5 0.2 \
--out ${qcdir}/${data}.LD

# Prune variants
plink2 \
--bfile ${qcdir}/${data}.IBD \
--allow-extra-chr \
--extract ${qcdir}/${data}.LD.prune.in \
--make-bed \
--out ${qcdir}/${data}.LD

# PCA on LD pruned, IBD and sib checked data
 
 plink2 \
 --bfile ${qcdir}/${data}.LD \
 --pca \
 --out ${qcdir}/${data}.LD 

mv *.eigenval ${pcadir}
mv *.eigenvec ${pcadir}
mv *.log ${log}